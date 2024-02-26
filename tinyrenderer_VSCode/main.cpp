#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <iostream>

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "my_gl.h"

Model *model = NULL;
int *zbuffer = NULL;
Vec3f lightDir = Vec3f(1, 1, 1).normalize();
Vec3f cameraPos(1, 1, 3);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);

const int width = 800;
const int height = 800;

// 高洛德着色
struct GouraudShader : public IShader
{
    Vec3f varying_intensity;

    virtual Vec4f vertex(int iface, int nthvert)
    {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert));                              // read the vertex from .obj file
        gl_Vertex = Viewport * Projection * ModelView * gl_Vertex;                            // transform it to screen coordinates
        varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert) * lightDir); // get diffuse lighting intensity
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color)
    {
        float intensity = varying_intensity * bar;
        if (intensity > .85)
            intensity = 1;
        else if (intensity > .60)
            intensity = .80;
        else if (intensity > .45)
            intensity = .60;
        else if (intensity > .30)
            intensity = .45;
        else if (intensity > .15)
            intensity = .30;
        else
            intensity = 0;
        color = TGAColor(255, 155, 0) * intensity;
        return false;
    }
};

// 漫反射
struct DiffuseShader : public IShader
{
    Vec3f varying_intensity;     // written by vertex shader, read by fragment shader
    mat<2, 3, float> varying_uv; // same as above

    virtual Vec4f vertex(int iface, int nthvert)
    {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert) * lightDir); // get diffuse lighting intensity
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert));                              // read the vertex from .obj file
        return Viewport * Projection * ModelView * gl_Vertex;                                 // transform it to screen coordinates
    }

    virtual bool fragment(Vec3f bar, TGAColor &color)
    {
        float intensity = varying_intensity * bar; // interpolate intensity for the current pixel
        Vec2f uv = varying_uv * bar;               // interpolate uv for the current pixel
        color = model->diffuse(uv) * intensity;    // well duh
        return false;                              // no, we do not discard this pixel
    }
};

// Phong氏着色
struct PhongShader : public IShader
{
    mat<2, 3, float> varying_uv; // same as above
    mat<4, 4, float> uniform_M = Projection * ModelView;
    mat<4, 4, float> uniform_MIT = ModelView.invert_transpose();
    virtual Vec4f vertex(int iface, int nthvert)
    {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        return Viewport * Projection * ModelView * gl_Vertex;    // transform it to screen coordinates
    }
    virtual bool fragment(Vec3f bar, TGAColor &color)
    {
        Vec2f uv = varying_uv * bar;
        Vec3f n = proj<3>(uniform_MIT * embed<4>(model->normal(uv))).normalize();
        Vec3f l = proj<3>(uniform_M * embed<4>(lightDir)).normalize();
        Vec3f r = (n * (n * l * 2.f) - l).normalize(); // reflected light
        float spec = pow(std::max(r.z, 0.0f), model->specular(uv));
        float diff = std::max(0.f, n * l);
        TGAColor c = model->diffuse(uv);
        color = c;
        for (int i = 0; i < 3; i++)
            color[i] = std::min<float>(5 + c[i] * (diff + .6 * spec), 255);
        return false;
    }
};

int main(int argc, char **argv)
{
    if (2 == argc)
    {
        model = new Model(argv[1]);
    }
    else
    {
        model = new Model("obj/african_head.obj");
    }

    lookat(cameraPos, center, up);
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    projection(-1.f / (cameraPos - center).norm());
    lightDir.normalize();

    TGAImage image(width, height, TGAImage::RGB);
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

    PhongShader shader;
    for (int i = 0; i < model->nfaces(); i++)
    {
        Vec4f screen_coords[3];
        for (int j = 0; j < 3; j++)
        {
            screen_coords[j] = shader.vertex(i, j);
        }
        triangle(screen_coords, shader, image, zbuffer);
    }

    image.flip_vertically(); // to place the origin in the bottom left corner of the image
    zbuffer.flip_vertically();
    image.write_tga_file("output.tga");
    zbuffer.write_tga_file("zbuffer.tga");

    delete model;
    return 0;
}
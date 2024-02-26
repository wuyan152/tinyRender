#ifndef __MY_GL_H__
#define __MY_GL_H__
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"

extern Matrix ModelView;
extern Matrix Viewport;
extern Matrix Projection;
const float depth = 255.f;

Vec3i m2v(Vec4f m);
Matrix v2m(Vec3f v);
void viewport(int x, int y, int w, int h);
void projection(float coeff = 0.f); // coeff = -1/c
void lookat(Vec3f eye, Vec3f center, Vec3f up);

struct IShader
{
    virtual ~IShader();
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};

// void triangle(Vec3i t[3], Vec2i uv[3], float intensity[3], TGAImage &image, int *zbuffer, int width, int height, Model *model);
void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer);
#endif //__MY_GL_H__
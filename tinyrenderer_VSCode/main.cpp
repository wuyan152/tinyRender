#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);

Model *model = NULL;
int *zbuffer = NULL;
Vec3f lightDir(0, 0, -1);
Vec3f cameraPos(0, 0, 3);

const int width = 800;
const int height = 800;
const int depth = 255;

// 叉乘
Vec3f crossProduct(const Vec3f &a, const Vec3f &b)
{
    Vec3f result;
    Vec3f _a(a), _b(b); // 创建新的 Vec3f 对象
    result[0] = _a[1] * _b[2] - _a[2] * _b[1];
    result[1] = _a[2] * _b[0] - _a[0] * _b[2];
    result[2] = _a[0] * _b[1] - _a[1] * _b[0];
    return result;
}

// 透视除法
Vec3f m2v(Matrix m)
{
    return Vec3f(m[0][0] / m[3][0], m[1][0] / m[3][0], m[2][0] / m[3][0]);
}

// 将三维坐标转换为4 x 1的矩阵
Matrix v2m(Vec3f v)
{
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

Matrix viewport(int x, int y, int w, int h)
{
    // 初始化一个4 x 4的单位矩阵
    Matrix m = Matrix::identity(4);
    // 以当前三维空间的中心为坐标原点
    m[0][3] = x + w / 2.f;
    m[1][3] = y + h / 2.f;
    m[2][3] = depth / 2.f;

    // 进行缩放操作，将[0,1]的坐标系缩放到[-0.5,0.5]，这样原点就被设置在了视口的中心
    m[0][0] = w / 2.f;
    m[1][1] = h / 2.f;
    m[2][2] = depth / 2.f;
    return m;
}

// //画线
// void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color)
// {
//     // 对于斜率超过1的线，我们交换x和y的坐标，使得之后使用x进行循环的时候，能打出最大数量的点，如果不这么做，y轴上可能会显得跳跃太大而产生孔洞
//     bool steep = false;
//     if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y))
//     {
//         std::swap(p0.x, p0.y);
//         std::swap(p1.x, p1.y);
//         steep = true;
//     }
//     // 一律从左往右绘制
//     if (p0.x > p1.x)
//     {
//         std::swap(p0, p1);
//     }
//     // 绘制循环
//     for (int x = p0.x; x <= p1.x; x++)
//     {
//         float t = (x - p0.x) / (float)(p1.x - p0.x);
//         int y = p0.y * (1. - t) + p1.y * t;
//         if (steep)
//         {
//             image.set(y, x, color);
//         }
//         else
//         {
//             image.set(x, y, color);
//         }
//     }
// }

// // 2d世界空间转屏幕空间
// Vec3f world2screen(Vec3f v)
// {
//     return Vec3f(int((v.x + 1.) * width / 2. + .5), int((v.y + 1.) * height / 2. + .5), v.z);
// }

// 求重心坐标
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P)
{
    Vec3f s[2];
    for (int i = 2; i--;)
    {
        s[i][0] = C[i] - A[i];
        s[i][1] = B[i] - A[i];
        s[i][2] = A[i] - P[i];
    }
    Vec3f u = crossProduct(s[0], s[1]);
    if (std::abs(u[2]) > 1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
    return Vec3f(-1, 1, 1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void triangle(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage &image, float intensity, int *zbuffer)
{
    if (t0.y == t1.y && t0.y == t2.y)
        return; // i dont care about degenerate triangles
    if (t0.y > t1.y)
    {
        std::swap(t0, t1);
        std::swap(uv0, uv1);
    }
    if (t0.y > t2.y)
    {
        std::swap(t0, t2);
        std::swap(uv0, uv2);
    }
    if (t1.y > t2.y)
    {
        std::swap(t1, t2);
        std::swap(uv1, uv2);
    }

    int total_height = t2.y - t0.y;
    for (int i = 0; i < total_height; i++)
    {
        bool second_half = i > t1.y - t0.y || t1.y == t0.y;
        int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
        float alpha = (float)i / total_height;
        float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here
        Vec3i A = t0 + Vec3f(t2 - t0) * alpha;
        Vec3i B = second_half ? t1 + Vec3f(t2 - t1) * beta : t0 + Vec3f(t1 - t0) * beta;
        Vec2i uvA = uv0 + (uv2 - uv0) * alpha;
        Vec2i uvB = second_half ? uv1 + (uv2 - uv1) * beta : uv0 + (uv1 - uv0) * beta;
        if (A.x > B.x)
        {
            std::swap(A, B);
            std::swap(uvA, uvB);
        }
        for (int j = A.x; j <= B.x; j++)
        {
            float phi = B.x == A.x ? 1. : (float)(j - A.x) / (float)(B.x - A.x);
            Vec3i P = Vec3f(A) + Vec3f(B - A) * phi;
            Vec2i uvP = uvA + (uvB - uvA) * phi;
            int idx = P.x + P.y * width;
            if (zbuffer[idx] < P.z)
            {
                zbuffer[idx] = P.z;
                TGAColor color = model->diffuse(uvP);
                image.set(P.x, P.y, TGAColor(color.r * intensity, color.g * intensity, color.b * intensity));
            }
        }
    }
}

// 使用重心坐标差值计算uv的方法，由于经过了投影变化后（非线性变换）后的坐标无法计算出正确的结果，故暂时弃用
// void triangle(Vec3f *pts, Vec2i *uv, int *zbuffer, TGAImage &image, float intensity)
// {
//     Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
//     Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
//     Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
//     for (int i = 0; i < 3; i++)
//     {
//         for (int j = 0; j < 2; j++)
//         {
//             bboxmin[j] = std::max(0.f, std::min(bboxmin[j], pts[i][j]));
//             bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
//         }
//     }
//     Vec3i P;
//     // std::cout << bboxmin.x << " " << bboxmin.y << " " << bboxmax.x << " " << bboxmax.y << " " << std::endl;
//     for (P.x = bboxmin[0]; P.x <= bboxmax[0]; P.x++)
//     {
//         for (P.y = bboxmin[1]; P.y <= bboxmax[1]; P.y++)
//         {
//             Vec3i bc_screen = barycentric(pts[0], pts[1], pts[2], P);
//             if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0)
//                 continue;
//             P.z = 0;
//             Vec2i uvP = Vec2i(0, 0);

//             P.z = pts[0].z * bc_screen.x + pts[1].z * bc_screen.y + pts[2].z * bc_screen.z;
//             uvP = uv[0] * bc_screen.x + uv[1] * bc_screen.y + uv[2] * bc_screen.z;

//             std::cout << P.z << " " << uvP.x << " " << uvP.y << std::endl;

//             int idx = int(P.x + P.y * width);
//             if (zbuffer[idx] < P.z)
//             {
//                 zbuffer[idx] = P.z;
//                 TGAColor color = model->diffuse(uvP);
//                 image.set(P.x, P.y, TGAColor(color.r * intensity, color.g * intensity, color.b * intensity, 255));
//             }
//         }
//     }
// }

int main(int argc, char **argv)
{

    // 如果不指定渲染的模型，则默认是else里的模型
    if (2 == argc)
    {
        model = new Model(argv[1]);
    }
    else
    {
        model = new Model("obj/african_head.obj");
    }

    zbuffer = new int[width * height];
    for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<int>::max())
        ;

    TGAImage image(width, height, TGAImage::RGB);

    // 透视投影矩阵
    Matrix Projection = Matrix::identity(4);
    // 视口变换矩阵
    Matrix ViewPort = viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    Projection[3][2] = -1.f / cameraPos.z;

    for (int i = 0; i < model->nfaces(); i++)
    {
        std::vector<int> face = model->face(i); // 获取模型的第i个面片
        Vec3f screenCoords[3];                  // 存贮第i个面片三个顶点的屏幕坐标
        Vec3f worldCoords[3];                   // 存储第I个面片三个顶点的世界坐标
        for (int j = 0; j < 3; j++)
        {
            Vec3f v = model->vert(face[j]);
            screenCoords[j] = m2v(ViewPort * Projection * v2m(v));
            worldCoords[j] = v;
        }

        Vec3f normal = crossProduct((worldCoords[2] - worldCoords[0]), (worldCoords[1] - worldCoords[0]));
        normal.normalize();
        float intensity = normal * lightDir;
        if (intensity > 0)
        {
            Vec2i uv[3];
            for (int j = 0; j < 3; j++)
                uv[j] = model->uv(i, j);
            // triangle(screenCoords, uv, zbuffer, image, intensity);
            triangle(screenCoords[0], screenCoords[1], screenCoords[2], uv[0], uv[1], uv[2], image, intensity, zbuffer);
        }
    }

    image.flip_vertically(); // 以左下角为坐标原点
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}
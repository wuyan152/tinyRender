#include <cmath>
#include <limits>
#include <cstdlib>
#include <iostream>

#include "my_gl.h"

Matrix ModelView;
Matrix Viewport;
Matrix Projection;

IShader::~IShader() {}

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
Vec3i m2v(Vec4f m)
{
    return Vec3i(m[0] / m[3], m[1] / m[3], m[2] / m[3]);
}

// 将三维坐标转换为4 x 4的单位矩阵
Matrix v2m(Vec3f v)
{
    Matrix m = Matrix::identity();
    m[0][3] = v.x;
    m[1][3] = v.y;
    m[2][3] = v.z;
    m[3][3] = 1.f;
    return m;
}

// 视口变换
void viewport(int x, int y, int w, int h)
{
    // 初始化一个4 x 4的单位矩阵
    Viewport = Matrix::identity();
    // 以当前三维空间的中心为坐标原点
    Viewport[0][3] = x + w / 2.f;
    Viewport[1][3] = y + h / 2.f;
    Viewport[2][3] = depth / 2.f;

    // 进行缩放操作，将[0,1]的坐标系缩放到[-0.5,0.5]，这样原点就被设置在了视口的中心
    Viewport[0][0] = w / 2.f;
    Viewport[1][1] = h / 2.f;
    Viewport[2][2] = depth / 2.f;
}

// 投影矩阵
void projection(float coeff)
{
    Projection = Matrix::identity();
    Projection[3][2] = coeff;
}

// 将世界空间坐标系转换到摄像机坐标系下
void lookat(Vec3f eye, Vec3f center, Vec3f up)
{
    Vec3f z = (eye - center).normalize();
    Vec3f x = crossProduct(up, z).normalize();
    Vec3f y = crossProduct(z, x).normalize();
    ModelView = Matrix::identity();
    for (int i = 0; i < 3; i++)
    {
        ModelView[0][i] = x[i];
        ModelView[1][i] = y[i];
        ModelView[2][i] = z[i];
        ModelView[i][3] = -center[i];
    }
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
Vec3f barycentric(Vec2f A, Vec2f B, Vec2f C, Vec2f P)
{
    Vec3f s[2];
    for (int i = 2; i--;)
    {
        s[i][0] = C[i] - A[i];
        s[i][1] = B[i] - A[i];
        s[i][2] = A[i] - P[i];
    }
    Vec3f u = cross(s[0], s[1]);
    if (std::abs(u[2]) > 1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
    return Vec3f(-1, 1, 1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
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

void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer)
{
    // 计算出三角包围盒
    Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            bboxmin[j] = std::min(bboxmin[j], pts[i][j] / pts[i][3]);
            bboxmax[j] = std::max(bboxmax[j], pts[i][j] / pts[i][3]);
        }
    }
    // 当前像素点
    Vec2i P;
    // 最终像素颜色
    TGAColor color;
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++)
    {
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++)
        {
            Vec3f c = barycentric(proj<2>(pts[0] / pts[0][3]), proj<2>(pts[1] / pts[1][3]), proj<2>(pts[2] / pts[2][3]), proj<2>(P));
            float z = pts[0][2] * c.x + pts[1][2] * c.y + pts[2][2] * c.z;
            float w = pts[0][3] * c.x + pts[1][3] * c.y + pts[2][3] * c.z;
            int frag_depth = std::max(0, std::min(255, int(z / w + .5)));
            if (c.x < 0 || c.y < 0 || c.z < 0 || zbuffer.get(P.x, P.y)[0] > frag_depth)
                continue;
            bool discard = shader.fragment(c, color);
            if (!discard)
            {
                zbuffer.set(P.x, P.y, TGAColor(frag_depth));
                image.set(P.x, P.y, color);
            }
        }
    }
}
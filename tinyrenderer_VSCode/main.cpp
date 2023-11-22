#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);

Model *model = NULL;
const int width = 800;
const int height = 800;

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color)
{
    // 对于斜率超过1的线，我们交换x和y的坐标，使得之后使用x进行循环的时候，能打出最大数量的点，如果不这么做，y轴上可能会显得跳跃太大而产生孔洞
    bool steep = false;
    if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y))
    {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }
    // 一律从左往右绘制
    if (p0.x > p1.x)
    {
        std::swap(p0, p1);
    }
    // 绘制循环
    for (int x = p0.x; x <= p1.x; x++)
    {
        float t = (x - p0.x) / (float)(p1.x - p0.x);
        int y = p0.y * (1. - t) + p1.y * t;
        if (steep)
        {
            image.set(y, x, color);
        }
        else
        {
            image.set(x, y, color);
        }
    }
}

// void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color)
// {
//     if (t0.y == t1.y && t0.y == t2.y)
//         return; // i dont care about degenerate triangles
//     if (t0.y > t1.y)
//         std::swap(t0, t1);
//     if (t0.y > t2.y)
//         std::swap(t0, t2);
//     if (t1.y > t2.y)
//         std::swap(t1, t2);
//     int total_height = t2.y - t0.y;
//     for (int i = 0; i < total_height; i++)
//     {
//         bool second_half = i > t1.y - t0.y || t1.y == t0.y;
//         int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
//         float alpha = (float)i / total_height;
//         float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here
//         Vec2i A = t0 + (t2 - t0) * alpha;
//         Vec2i B = second_half ? t1 + (t2 - t1) * beta : t0 + (t1 - t0) * beta;
//         if (A.x > B.x)
//             std::swap(A, B);
//         for (int j = A.x; j <= B.x; j++)
//         {
//             image.set(j, t0.y + i, color); // attention, due to int casts t0.y+i != A.y
//         }
//     }
// }

// Vec3f barycentric(Vec2i *pts, Vec2i P) {
//     Vec3f u = Vec3f(pts[2][0]-pts[0][0], pts[1][0]-pts[0][0], pts[0][0]-P[0])^Vec3f(pts[2][1]-pts[0][1], pts[1][1]-pts[0][1], pts[0][1]-P[1]);
//     /* `pts` and `P` has integer value as coordinates
//        so `abs(u[2])` < 1 means `u[2]` is 0, that means
//        triangle is degenerate, in this case return something with negative coordinates */
//     if (std::abs(u.z)<1) return Vec3f(-1,1,1);
//     return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
// }

// void triangle(Vec2i *pts, TGAImage &image, TGAColor color) {
//     Vec2i bboxmin(image.get_width()-1,  image.get_height()-1);
//     Vec2i bboxmax(0, 0);
//     Vec2i clamp(image.get_width()-1, image.get_height()-1);
//     for (int i=0; i<3; i++) {
//         bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
// 	bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));

// 	bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
// 	bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
//     }
//     Vec2i P;
//     for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
//         for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
//             Vec3f bc_screen  = barycentric(pts, P);
//             if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
//             image.set(P.x, P.y, color);
//         }
//     }
// }

/**三角形的光栅化，个人思考*/
// 二维叉乘
float crossProduct(Vec2i a, Vec2i b)
{
    return a.x * b.y - a.y * b.x;
}

// 判断一个点是否在三角形内部
bool isInTriangle(Vec2i t0, Vec2i t1, Vec2i t2, Vec2i point)
{
    Vec2i v0 = t1 - t0;
    Vec2i v1 = t2 - t1;
    Vec2i v2 = t0 - t2;

    Vec2i v0P = point - t0;
    Vec2i v1P = point - t1;
    Vec2i v2P = point - t2;

    bool t0P = crossProduct(v0, v0P) >= 0;
    bool t1P = crossProduct(v1, v1P) >= 0;
    bool t2P = crossProduct(v2, v2P) >= 0;

    return t0P == t1P && t1P == t2P && t0P == t2P;
}

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color)
{
    int tranMaxX = std::max(t0.x, std::max(t1.x, t2.x));
    int tranMinX = std::min(t0.x, std::min(t1.x, t2.x));
    int tranMaxY = std::max(t0.y, std::max(t1.y, t2.y));
    int tranMinY = std::min(t0.y, std::min(t1.y, t2.y));

    for (int i = tranMinX; i <= tranMaxX; i++)
    {
        for (int j = tranMinY; j <= tranMaxY; j++)
        {
            if (isInTriangle(t0, t1, t2, Vec2i(i, j)))
            {
                image.set(i, j, color);
            }
        }
    }
}

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

    TGAImage image(width, height, TGAImage::RGB);

    Vec3f light_dir(0, 0, -1);
    for (int i = 0; i < model->nfaces(); i++)
    {
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3];
        Vec3f world_coords[3];
        for (int j = 0; j < 3; j++)
        {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
            world_coords[j] = v;
        }
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        n.normalize();
        float intensity = n * light_dir;
        if (intensity > 0)
        {
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
        }
    }
    image.flip_vertically(); // 以左下角为坐标原点
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}

#/** * * * BME-VIK-MI_2018_02_08 * * *|3.HF| * * * * * * * * * *\
#*    _ _____   _____        __ _                              *
#*   (_)  __ \ / ____|      / _| |                             *
#*   |_| |__)| (___    ___ | |_| |___      ____ _ _ __ ___     *
#*   | |  _  / \___ \ / _ \|  _| __\ \ /\ / / _` | '__/ _ \    *
#*   | | | \ \ ____) | (_) | | | |_ \ V  V / (_| | | |  __/    *
#*   |_|_|  \_\_____/ \___/|_|  \__| \_/\_/ \__,_|_|  \___|    *
#*                                                             *
#*                   http://irsoftware.net                     *
#*                                                             *
#*              contact_adress: sk8Geri@gmail.com               *
#*                                                               *
#*       This file is a part of the work done by aFagylaltos.     *
#*         You are free to use the code in any way you like,      *
#*         modified, unmodified or copied into your own work.     *
#*        However, I would like you to consider the following:    *
#*                                                               *
#*  -If you use this file and its contents unmodified,         *
#*              or use a major part of this file,               *
#*     please credit the author and leave this note untouched.   *
#*  -If you want to use anything in this file commercially,      *
#*                please request my approval.                    *
#*                                                              *
#\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include <mutex>

#if defined(__APPLE__)

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#else
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

const int screenWidth = 600;
const int screenHeight = 600;
const float E = 0.001f;
float image[screenWidth * screenHeight * 3];
const size_t numBubies = 4;
std::mutex mutex;
std::atomic<bool> running;

enum CSG_operator {
    UNIO,
    FELTER_TOP,
    FELTER_BOTTOM
};

float randFloat(float a, float b) {
    return ((b - a) * ((float) rand() / RAND_MAX)) + a;
}

//--------------------------------------------------------
// 3D Vektor
//--------------------------------------------------------
struct Vector {
    float x, y, z, h;

    Vector() {
        x = y = z = h = 0.0f;
    }

    Vector(float x0, float y0, float z0, float h0 = 0) {
        x = x0;
        y = y0;
        z = z0;
        h = h0;
    }

    float &operator[](int i) { return *(&x + i); }

    Vector operator*(float a) const {
        return Vector(x * a, y * a, z * a, h * a);
    }

    Vector operator/(float d) const {
        return Vector(x / d, y / d, z / d, h / d);
    }

    Vector operator+(const Vector &v) const {
        return Vector(x + v.x, y + v.y, z + v.z, h + v.h);
    }

    Vector operator-(const Vector &v) const {
        return Vector(x - v.x, y - v.y, z - v.z, h - v.h);
    }

    float operator*(const Vector &v) const {    // dot product
        return (x * v.x + y * v.y + z * v.z + h * v.h);
    }

    Vector operator%(const Vector &v) {
        return Vector(y * v.z - v.y * z, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    float Length() { return sqrtf(x * x + y * y + z * z + h * h); }

    Vector &Normalize() {
        float norm = Length();
        x /= norm;
        y /= norm;
        z /= norm;
        h /= norm;
        return *this;
    }
};

#define Point(x, y, z) Vector(x,y,z,1)

struct Matrix4x4 {
    float m[4][4];

    Matrix4x4() {
        clear();
    }

    void clear() {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                m[i][j] = 0;
    }

    Vector operator*(Vector v) {
        Vector result(0.0f, 0.0f, 0.0f, 0.0f);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result[i] += m[i][j] * v[j];

        return result;
    }

    Matrix4x4 operator*(const Matrix4x4 &mat) {
        Matrix4x4 result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                result.m[i][j] = 0;
                for (int k = 0; k < 4; k++)
                    result.m[i][j] += m[i][k] * mat.m[k][j];
            }

        return result;
    }

    float &operator()(int i, int j) { return m[i][j]; }
};

struct Color {
    float r, g, b;

    Color() {
        r = g = b = 0;
    }

    Color(float r0, float g0, float b0) {
        r = r0;
        g = g0;
        b = b0;
    }

    Color operator*(float a) {
        return Color(r * a, g * a, b * a);
    }

    Color operator*(const Color &c) {
        return Color(r * c.r, g * c.g, b * c.b);
    }

    Color operator+(const Color &c) {
        return Color(r + c.r, g + c.g, b + c.b);
    }

    Color operator-(const Color &c) {
        return Color(r - c.r, g - c.g, b - c.b);
    }

    Color operator+=(const Color &c) {
        r += c.r;
        g += c.g;
        b += c.b;
        return *this;
    }
};

struct Ray {
    Vector o;
    Vector v;
};

class Camera {
public:
    Vector eye;
    Vector lookat;
    Vector up;
    Vector right;
public:
    Camera(Vector _eye, Vector _lookat, Vector _right, Vector _up) {
        eye = _eye;
        lookat = _lookat;
        right = _right;
        up = _up;
    }

    Ray getRay(int x, int y) {
        Ray ray;
        ray.o = eye;
        ray.v = lookat + right * (2.0f * x / screenWidth - 1) + up * (2.0f * y / screenHeight - 1) - eye;
        ray.v.Normalize();
        return ray;
    }
};

struct Material {
    Color F0;
    Color ka, kd, ks;
    float n, shine;
    bool isReflective, isRefractive;

    Material() {
        isReflective = false;
        isRefractive = false;
    }

    void ReflectionDir(Vector &R, Vector &N, Vector &V) {
        float cosa = -(N * V);
        R = V + N * cosa * 2;
    }

    bool RefractionDir(Vector &T, Vector &N, Vector &V) {
        float cosa = -(N * V), cn = n;
        if (cosa < 0) {
            cosa = -cosa;
            N = N * (-1);
            cn = 1 / n;
        }
        float disc = 1 - (1 - cosa * cosa) / cn / cn;
        if (disc < 0) return false;
        T = V / cn + N * (cosa / cn - sqrt(disc));
        return true;
    }

    Color Frensel(Vector &N, Vector &V) {
        float cosa = fabs(N * V);
        return F0 + (Color(1, 1, 1) - F0) * pow(1 - cosa, 5);
    }
};

struct Hit {
    Vector x;
    Vector normal;
    Material mat;
    float t;

    Hit() {
        t = -1;
        actual = false;
    }

    bool actual;

    void select() {
        actual = true;
    }

    void deselect() {
        actual = false;
    }

    bool isActual() {
        return actual;
    }
};

class Object {
protected:
    Material mat;
public:
    Object(Material _mat) {
        mat = _mat;
    }

    Material getMaterial() {
        return mat;
    }

    int indexOfSelected(Hit *hits[], int numHit) {
        for (int i = 0; i < numHit; i++) {
            if (hits[i]->isActual())
                return i;
        }
        return -1;
    }

    int hitShort(Hit *hit, Hit *hits[], int numHit, int csg_operator) {
        hit->select();

        if (numHit == 0)
            hits[numHit++] = hit;
        else if (hits[0] != NULL && hits[numHit - 1]->t < hit->t) {
            hits[numHit++] = hit;
        } else if (hits[0] != NULL) {
            for (int i = 0; i < numHit; i++) {
                if (hits[i]->t > hit->t) {

                    for (int j = numHit; j > i; j--)
                        hits[j] = hits[j - 1];

                    hits[i] = hit;
                    numHit++;
                    break;
                }
            }
        }

        if (hits[0] == NULL)
            return 0;

        int selected = indexOfSelected(hits, numHit);

        if (csg_operator == FELTER_TOP) {

            if (numHit > 1 && selected == 0) {
                delete hits[selected];
                hits[selected] = NULL;
            }

            if (hits[selected + 1] == NULL) {
                delete hits[selected];
                hits[selected] = NULL;
            }

            for (int i = 0; i < numHit; i++) {
                if (hits[i] != NULL)
                    if (hits[i]->x.x > 4.5 && i != selected) {
                        delete hits[i];
                        hits[i] = NULL;
                    }
            }
        } else if (csg_operator == FELTER_BOTTOM) {
            for (int i = 0; i < numHit; i++) {
                if (hits[i]->x.x < -0.01) {
                    delete hits[i];
                    hits[i] = NULL;
                }
            }
        }

        for (int i = 0; i < numHit; i++) {
            if (hits[i] == NULL) {

                for (int j = i; j < numHit; j++) {
                    hits[j] = hits[j + 1];
                    hits[j + 1] = NULL;
                }

                numHit--;
                i--;
            }
        }

        if (selected != -1 && hits[selected] != NULL)
            hits[selected]->deselect();

        return numHit;
    }

    virtual Hit Intersect(Ray ray) = 0;

    virtual void CSGIntersect(Ray ray, Hit *hits[], int &numHit, int csg_operator) = 0;
};

class Elipszoid : public Object {
    Matrix4x4 A;
public:
    Elipszoid(Vector position, Vector axes, Material mat) : Object(mat) {
        float a = axes[0];
        float b = axes[1];
        float c = axes[2];

        float a2 = a * a;
        float b2 = b * b;
        float c2 = c * c;

        float x0 = position[0];
        float y0 = position[1];
        float z0 = position[2];

        A(0, 0) = 1.0f / a2;
        A(0, 1) = 0;
        A(0, 2) = 0;
        A(0, 3) = -x0 / a2;
        A(1, 0) = 0;
        A(1, 1) = 1.0f / b2;
        A(1, 2) = 0;
        A(1, 3) = -y0 / b2;
        A(2, 0) = 0;
        A(2, 1) = 0;
        A(2, 2) = 1.0f / c2;
        A(2, 3) = -z0 / c2;
        A(3, 0) = -x0 / a2;
        A(3, 1) = -y0 / b2;
        A(3, 2) = -z0 / c2;
        A(3, 3) = -1 + x0 * x0 / a2 + y0 * y0 / b2 + z0 * z0 / c2;

    }

    Hit Intersect(Ray ray) override {
        Hit hit;
        double a = ray.v * (A * ray.v);
        double b = ray.o * (A * ray.v) + ray.v * (A * ray.o);
        double c = ray.o * (A * ray.o);
        double disc = b * b - 4 * a * c;
        if (disc < E) return hit;
        disc = sqrt(disc);
        hit.t = (-b - disc) / 2 / a;

        if (hit.t < E) hit.t = (-b + disc) / 2 / a;
        if (hit.t < E) return hit;
        hit.x = ray.o + ray.v * hit.t;
        hit.normal = A * hit.x;
        hit.normal.h = 0;
        hit.mat = mat;
        return hit;
    }

    void CSGIntersect(Ray ray, Hit *hits[], int &numHit, int csg_operator) override {
        double a = ray.v * (A * ray.v);
        double b = ray.o * (A * ray.v) + ray.v * (A * ray.o);
        double c = ray.o * (A * ray.o);

        double disc = b * b - 4 * a * c;

        if (disc < E) {
            return;
        }

        disc = sqrt(disc);
        Hit *hit = new Hit();
        hit->deselect();
        hit->t = (-b - disc) / 2 / a;
        if (hit->t >= E) {
            hit->x = ray.o + ray.v * hit->t;
            hit->normal = A * hit->x;
            hit->normal.h = 0;
            hit->mat = mat;

            numHit = hitShort(hit, hits, numHit, csg_operator);
        } else delete hit;

        hit = new Hit();
        hit->deselect();
        hit->t = (-b + disc) / 2 / a;
        if (hit->t >= E) {
            hit->x = ray.o + ray.v * hit->t;
            hit->normal = A * hit->x;
            hit->normal.h = 0;
            hit->mat = mat;

            numHit = hitShort(hit, hits, numHit, csg_operator);
        } else delete hit;
    }
};

class Plane : public Object {
    Vector p1, p2, p3, p4;
    Vector n;
    Material mat1;
    Material mat2;

public:
    Plane(Vector _p1, Vector _p2, Vector _p3, Vector _p4, Material material1, Material material2)
            : Object(material1) {
        p1 = _p1;
        p2 = _p2;
        p3 = _p3;
        p4 = _p4;

        n = (p1 - p2) % (p3 - p2);
        n.Normalize();
        n.h = 0;

        mat1 = material1;
        mat2 = material2;
    }

    Hit Intersect(Ray ray) override {
        Hit hit;

        float a = ray.v * n;
        if (a < E && a > -E)
            return hit;

        hit.t = (p2 - ray.o) * n / a;

        if (hit.t < E)
            return hit;

        hit.x = ray.o + ray.v * hit.t;

        if (ray.o.x > p1.x)
            hit.normal = n * -1;
        else
            hit.normal = n;

        hit.mat = mat1;

        if ((((p1 - p2) % ((hit).x - p2)) * n >= 0 && ((p3 - p1) % ((hit).x - p1)) * n >= 0 &&
             ((p2 - p3) % ((hit).x - p3)) * n >= 0) ||
            (((p1 - p3) % ((hit).x - p3)) * n >= 0 && ((p4 - p1) % ((hit).x - p1)) * n >= 0 &&
             ((p3 - p4) % ((hit).x - p4)) * n >= 0)) {
            return hit;
        }

        hit.t = -1;
        return hit;
    }

    void CSGIntersect(Ray ray, Hit *hits[], int &numHit, int csg_operator) override {
        float a = ray.v * n;
        if (a < E && a > -E)
            return;

        Hit *hit = new Hit();
        hit->deselect();
        hit->t = (p2 - ray.o) * n / a;

        if (hit->t < E) {
            delete hit;
            return;
        }

        hit->x = ray.o + ray.v * hit->t;

        if (ray.o.x > p1.x)
            hit->normal = n * -1;
        else
            hit->normal = n;

        float yy = (hit->x.y);
        float zz = (hit->x.z);

        if (yy < 0) {
            if (((int) yy + (int) zz) % 2 == 0)
                hit->mat = mat1;
            else
                hit->mat = mat2;
        } else {
            if (((int) yy + (int) zz) % 2 == 0)
                hit->mat = mat2;
            else
                hit->mat = mat1;
        }

        Vector p = hit->x;

        if ((((p1 - p2) % (p - p2)) * n >= 0 && ((p3 - p1) % (p - p1)) * n >= 0 && ((p2 - p3) % (p - p3)) * n >= 0) ||
            (((p1 - p3) % (p - p3)) * n >= 0 && ((p4 - p1) % (p - p1)) * n >= 0 && ((p3 - p4) % (p - p4)) * n >= 0)) {
            numHit = hitShort(hit, hits, numHit, csg_operator);
            return;
        }

        delete hit;
        return;
    }
};

struct Light {
    Vector pos;
    Color lin;

    Light(Vector _pos, Color _lin) : pos(_pos), lin(_lin) {}
};

class CSG : public Object {
    Object *csg_node[8];
    int numNode;

    int csg_operator[9];
    int numOperator;
    int numHit;

public:
    CSG() : Object(Material()) {
        numNode = 0;
        numOperator = 0;
    }

    void addObject(Object *obj, int _csg_operator) {
        csg_node[numNode++] = obj;
        csg_operator[numOperator++] = _csg_operator;
    }

    void CSGIntersect(Ray ray, Hit *hits[], int &numHit, int csg_operator) override {}

    Hit Intersect(Ray ray) override {
        Hit *hits[10];
        numHit = 0;
        Hit hit;

        for (int i = 0; i < 10; i++) {
            hits[i] = NULL;
        }

        for (int i = 0; i < numNode; i++) {
            csg_node[i]->CSGIntersect(ray, hits, numHit, csg_operator[i]);
        }

        for (int i = 0; i < numHit; i++) {
            if (hits[i] != NULL) {
                hit = *hits[i];
                break;
            }
        }

        for (int i = 0; i < 10; i++) {
            if (hits[i] != 0)
                delete hits[i];
        }

        return hit;
    }
};

struct Sector {
    int x, y;
    int width, height;
    std::atomic<bool> finished;
    std::atomic<bool> rendering;

    Sector() : x(0), y(0), width(0), height(0), finished(false), rendering(false) {};

} sectors[8][8];

void sectorIJ(int &_i, int &_j) {
    std::lock_guard<std::mutex> guard(mutex);
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            if (!sectors[i][j].rendering) {
                sectors[i][j].rendering = true;
                _i = i;
                _j = j;
                return;
            }
        }
    }

    running.store(false);
}

struct Raytracer {
    Raytracer() {
        La = Color(0.2, 0.2, 0.5);
    }

    Color trace(const std::vector<Light> &lights, const std::vector<Object *> &objects, Ray ray, int d = 0);

    Hit intersectAll(const std::vector<Object *> &objects, Ray ray);

    Color La;
};

class Scene {
public:
    Camera camera;
    std::vector<Light> lights;
    std::vector<Object *> objects;
public:
    Scene();
    void makeDrawable();
    void Render();
};


Scene::Scene()
        : camera(Point(8, 0, -2), Point(7, 0, 0), Vector(0, 1, 0), Vector(1, 0, 0)) {

    lights.push_back(Light(Vector(2.5, 1.2, -0.2), Color(1.9, 1.9, 1.9)));
    lights.push_back(Light(Vector(2.5, 1.2, -0.2), Color(1.9, 1.8, 0.6)));
}

void Scene::makeDrawable() {
    objects.clear();

    Material glass_mat;
    glass_mat.kd = Color(0.01, 0.01, 0.02);
    glass_mat.ka = Color(0.5, 0.5, 0.5);
    glass_mat.ks = Color(0.0, 0.0, 0.0);
    glass_mat.shine = 50;
    glass_mat.n = 1.5f;
    glass_mat.F0 = Color(0.2, 0.2, 0.0);
    glass_mat.isReflective = false;
    glass_mat.isRefractive = true;

    Material bier_mat;
    bier_mat.kd = Color(0.0, 0.0, 0.0);
    bier_mat.ka = Color(0.3, 0.3, 0.0);
    bier_mat.ks = Color(0.0, 0.0, 0.0);
    bier_mat.shine = 50;
    bier_mat.n = 1.3 / 1.5;
    bier_mat.F0 = Color(0, 0, 0.8);
    bier_mat.isReflective = false;
    bier_mat.isRefractive = true;

    Material bier_mat2;
    bier_mat2.kd = Color(0.0, 0.0, 0.0);
    bier_mat2.ka = Color(0.4, 0.4, 0.0);
    bier_mat2.ks = Color(0.0, 0.0, 0.0);
    bier_mat2.shine = 50;
    bier_mat2.n = 1.3 / 1.0;
    bier_mat2.F0 = Color(0, 0, 0.9);
    bier_mat2.isReflective = true;
    bier_mat2.isRefractive = true;

    Material talp_mat;
    talp_mat.kd = Color(.1, .1, .1);
    talp_mat.ka = Color(.1, .1, .1);
    talp_mat.ks = Color(0.0, 0.0, 0.0);
    talp_mat.shine = 40;
    talp_mat.n = 1.5f;
    talp_mat.F0 = Color(0.0, 0.0, 0.0);
    talp_mat.isReflective = false;
    talp_mat.isRefractive = true;

    Material bubi_mat;
    bubi_mat.kd = Color(0.2, 0.5, 0.2);
    bubi_mat.ka = Color(0.3, 0.3, 0.3);
    bubi_mat.ks = Color(0.0, 0.0, 0.0);
    bubi_mat.shine = 20;
    bubi_mat.n = 1.0 / 1.3f;
    bubi_mat.F0 = Color(0.0, 0.0, 0.0);
    bubi_mat.isReflective = false;
    bubi_mat.isRefractive = true;

    Material plane_mat1;
    plane_mat1.kd = Color(1.0, 1.0, 1.0);
    plane_mat1.ka = Color(0.5, 0.5, 0.5);
    plane_mat1.ks = Color(4.0, 4.0, 4.0);
    plane_mat1.shine = 30;
    plane_mat1.n = 1.5f;
    plane_mat1.F0 = Color(0.5, 0.5, 0.5);
    plane_mat1.isReflective = false;
    plane_mat1.isRefractive = false;

    Material plane_mat2;
    plane_mat2.kd = Color(1.0, 0.0, 0.0);
    plane_mat2.ka = Color(0.8, 0.0, 0.0);
    plane_mat2.ks = Color(2.0, 2.0, 2.0);
    plane_mat2.shine = 30;
    plane_mat2.n = 1.5f;
    plane_mat2.F0 = Color(0.5, 0.5, 0.0);
    plane_mat2.isReflective = false;
    plane_mat2.isRefractive = false;

    Elipszoid *glass = new Elipszoid(Point(3.5, 0, 10), Vector(3, 1.5, 1.5), glass_mat);
    Elipszoid *bier = new Elipszoid(Point(3.5, 0, 10), Vector(2.9, 1.45, 1.45), bier_mat);
    Elipszoid *talp = new Elipszoid(Point(0, 0, 10), Vector(0.5, 1.2, 1.3), talp_mat);

    Plane *asztal = new Plane(Point(0, 10, 5), Point(0, 10, 25), Point(0, -10, 25), Point(0, -10, 5), plane_mat1,
                              plane_mat2);
    Plane *cutter = new Plane(Point(4.5, 10, 5), Point(4.5, 10, 30), Point(4.5, -10, 30), Point(4.5, -10, 5), bier_mat2,
                              bier_mat2);

    CSG *csg = new CSG();
    csg->addObject(bier, UNIO);
    csg->addObject(glass, UNIO);
    csg->addObject(cutter, FELTER_TOP);
    csg->addObject(talp, UNIO);
    csg->addObject(asztal, FELTER_BOTTOM);
    objects.push_back(csg);

    for (size_t i = 0; i < numBubies; i++) {
        float x = randFloat(-1.0, 1.0);
        float y = randFloat(1.0, 3.6);
        float z = randFloat(10.5, 11.0);

        float sizeX = randFloat(0.02, 0.05);
        float sizeY = randFloat(0.02, 0.04);
        float sizeZ = randFloat(0.02, 0.04);

        objects.push_back(new Elipszoid(Point(y, x, z), Vector(sizeX, sizeY, sizeZ), bubi_mat));
    }
}

void reset() {
    int w = 75;
    int h = 75;

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            sectors[i][j].x = i * w;
            sectors[i][j].y = j * h;
            sectors[i][j].width = w;
            sectors[i][j].height = h;
            sectors[i][j].finished = false;
            sectors[i][j].rendering = false;
        }
    }
}

void Scene::Render() {
    Raytracer raytracer;
    while(running)
    {
        int indexI = 0;
        int indexJ = 0;
        sectorIJ(indexJ, indexI);
        for (int x = sectors[indexI][indexJ].x; x < sectors[indexI][indexJ].width + sectors[indexI][indexJ].x; x++) {
            for (int y = sectors[indexI][indexJ].y;
                 y < sectors[indexI][indexJ].height + sectors[indexI][indexJ].y; y++) {
                Ray ray = camera.getRay(x, y);
                Color Lrad = raytracer.trace(lights, objects, ray);
                image[y * screenWidth * 3 + x * 3 + 0] = Lrad.r;
                image[y * screenWidth * 3 + x * 3 + 1] = Lrad.g;
                image[y * screenWidth * 3 + x * 3 + 2] = Lrad.b;
            }
        }
        std::lock_guard<std::mutex> guard(mutex);
        sectors[indexI][indexJ].finished = true;
    }
}

Color Raytracer::trace(const std::vector<Light> &lights, const std::vector<Object *> &objects, Ray ray, int d) {
    if (d > 5) return La;
    Hit hit = intersectAll(objects, ray);
    if (hit.t < 0) return La;

    Color c = La * hit.mat.ka;
    Vector N = hit.normal;
    Vector x = hit.x;

    Color temp = hit.mat.F0;

    for (int l = 0; l < lights.size(); l++) {
        Ray shadowRay;
        shadowRay.o = x;
        shadowRay.v = lights[l].pos - x;
        Hit shadowHit = intersectAll(objects, shadowRay);
        Vector y = shadowHit.x;
        if (shadowHit.t<0 || (x - y).Length()>(x - lights[l].pos).Length()) {
            Vector H = lights[l].pos - ray.v;
            H.Normalize();
            float costheta = (N * shadowRay.v.Normalize());
            if (costheta < 0) costheta = 0;
            float cosdelta = H * N;
            if (cosdelta < 0 /*|| hit.mat->isRefractive*/) cosdelta = 0;
            Color lin = lights[l].lin;
            c += lin * (hit.mat.kd * costheta + hit.mat.ks * pow(cosdelta, hit.mat.shine));
        }
    }

    if (hit.mat.isReflective) {
        Ray reflectedRay;
        reflectedRay.o = x;
        hit.mat.ReflectionDir(reflectedRay.v, N, ray.v);
        c += hit.mat.Frensel(N, ray.v) * trace(lights, objects, reflectedRay, d + 1);
    }

    if (hit.mat.isRefractive) {
        Ray refractedRay;
        refractedRay.o = x;

        if (hit.mat.RefractionDir(refractedRay.v, N, ray.v))
            c += (Color(1, 1, 1) - hit.mat.Frensel(N, ray.v)) * trace(lights, objects, refractedRay, d + 1);
        else {
            hit.mat.F0 = Color(1, 1, 1);
            Ray reflectedRay;
            reflectedRay.o = x;
            hit.mat.ReflectionDir(reflectedRay.v, N, ray.v);
            c += hit.mat.Frensel(N, ray.v) * trace(lights, objects, reflectedRay, d + 1);
        }
    }

    hit.mat.F0 = temp;

    return c;
}

Hit Raytracer::intersectAll(const std::vector<Object *> &objects, Ray ray) {
    Hit hit;
    for (int i = 0; i < objects.size(); i++) {
        Hit newHit = objects[i]->Intersect(ray);
        if (newHit.t > E) {
            if (hit.t < E || newHit.t < hit.t)
                hit = newHit;
        }
    }

    if (hit.t > E) hit.normal.Normalize();
    return hit;
}

#define THREADS 16
std::thread threads[THREADS];
Scene scene[THREADS];

void onInitialization() {
    running.store(true);
    glViewport(0, 0, screenWidth, screenHeight);
    reset();
    for (int i = 0; i < THREADS; ++i) {
        scene[i].makeDrawable();
    }
}

void onDisplay() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    for (int i = 0; i < THREADS; ++i) {
        threads[i] = std::thread([i]() {
            scene[i].lights[0].pos.y += 0.01f;
            scene[i].lights[0].pos.z += 0.01f;
            scene[i].Render();
        });
    }

    for (int i = 0; i < THREADS; ++i) {
        threads[i].join();
    }

    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, image);

    glutSwapBuffers();

    bool finished = true;
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            if (!sectors[i][j].finished) {
                finished = false;
                break;
            }
        }
    }

    if (finished) {
        std::cout << "finished" << std::endl;
        reset();
        running.store(true);
    }
}

void onKeyboard(unsigned char key, int x, int y) {
}

void onMouse(int button, int state, int x, int y) {
}

void onIdle() {
    glutPostRedisplay();
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Raytrace - lesson 04");

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    onInitialization();

    glutDisplayFunc(onDisplay);
    glutMouseFunc(onMouse);
    glutIdleFunc(onIdle);
    glutKeyboardFunc(onKeyboard);
    glutMainLoop();

    return 0;
} 

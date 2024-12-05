#include <iostream> 
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(double t) const { return Vec3(x / t, y / t, z / t); }
    Vec3 operator*(const Vec3& v) const { return Vec3(x * v.x, y * v.y, z * v.z); }

    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3 normalize() const {
        double len = std::sqrt(x * x + y * y + z * z);
        return *this / len;
    }
};

struct Ray {
    Vec3 origin, direction;

    Ray(const Vec3& origin, const Vec3& direction)
        : origin(origin), direction(direction.normalize()) {}
};

struct Sphere {
    Vec3 center;
    double radius;
    Vec3 color;    // Diffuse color
    double reflectivity; // Reflectivity factor

    Sphere(Vec3 c, double r, Vec3 clr, double refl)
        : center(c), radius(r), color(clr), reflectivity(refl) {}
};

bool hit_sphere(const Sphere& sphere, const Ray& ray, double& t) {
    Vec3 oc = ray.origin - sphere.center;
    double a = ray.direction.dot(ray.direction);
    double b = 2.0 * oc.dot(ray.direction);
    double c = oc.dot(oc) - sphere.radius * sphere.radius;
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return false;
    t = (-b - std::sqrt(discriminant)) / (2.0 * a);
    return true;
}

// Define a plane
struct Plane {
    Vec3 point; // A point on the plane
    Vec3 normal; // Normal vector to the plane

    Plane(Vec3 p, Vec3 n) : point(p), normal(n) {}
};

// Check for intersection between ray and plane
bool hit_plane(const Plane& plane, const Ray& ray, double& t) {
    double denominator = ray.direction.dot(plane.normal);
    if (std::fabs(denominator) > 1e-6) { // Avoid division by zero
        Vec3 p0l0 = plane.point - ray.origin;
        t = p0l0.dot(plane.normal) / denominator;
        return (t >= 0); // Check if the intersection is in front of the ray
    }
    return false;
}

Vec3 ray_color(const Ray& ray, const Vec3& light_position, const std::vector<Sphere>& spheres, const Plane& plane) {
    Vec3 light_color(1, 1, 1); // White light
    double closest_t = INFINITY;
    const Sphere* hit_sphere_ptr = nullptr;

    // Find the closest sphere that the ray hits
    for (const auto& sphere : spheres) {
        double t;
        if (hit_sphere(sphere, ray, t) && t < closest_t) {
            closest_t = t;
            hit_sphere_ptr = &sphere;
        }
    }

    // Check for plane intersection
    double plane_t;
    bool hit_plane_flag = hit_plane(plane, ray, plane_t);
    if (hit_plane_flag && plane_t < closest_t) {
        // If the ray hits the plane first, return the plane color (light gray)
        return Vec3(0.8, 0.8, 0.8); // Light gray plane color
    }

    if (hit_sphere_ptr) {
        const Sphere& sphere = *hit_sphere_ptr;
        Vec3 hit_point = ray.origin + ray.direction * closest_t;
        Vec3 normal = (hit_point - sphere.center).normalize();
        Vec3 light_dir = (light_position - hit_point).normalize();
        double diffuse = std::max(0.0, normal.dot(light_dir));
        return sphere.color * light_color * diffuse;
    }

    return Vec3(0, 0, 0); // Black background
}

void render_frame(const std::string& filename, int image_width, int image_height, const Vec3& light_position, const std::vector<Sphere>& spheres, const Plane& plane) {
    std::ofstream file(filename);
    file << "P3\n" << image_width << " " << image_height << "\n255\n";

    Vec3 origin(0, 0, 0);
    double viewport_height = 2.0;
    double viewport_width = 2.0 * image_width / image_height;
    double focal_length = 1.0;

    Vec3 horizontal(viewport_width, 0, 0);
    Vec3 vertical(0, viewport_height, 0);
    Vec3 lower_left_corner = origin - horizontal / 2 - vertical / 2 - Vec3(0, 0, focal_length);

    for (int j = image_height - 1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            double u = double(i) / (image_width - 1);
            double v = double(j) / (image_height - 1);
            Ray ray(origin, lower_left_corner + horizontal * u + vertical * v - origin);
            Vec3 color = ray_color(ray, light_position, spheres, plane);
            file << static_cast<int>(255.999 * color.x) << " "
                << static_cast<int>(255.999 * color.y) << " "
                << static_cast<int>(255.999 * color.z) << "\n";
        }
    }

    file.close();
}

int main() {
    const int image_width = 400;
    const int image_height = 400;
    const int frame_count = 60;
    const double radius = 2.0;

    // Define spheres with different materials
    std::vector<Sphere> spheres = {
        Sphere(Vec3(0, 0, -3), 0.5, Vec3(1, 1, 0), 0.5), // Yellow sphere
        Sphere(Vec3(-1, 0, -2), 0.5, Vec3(0, 1, 1), 0.3), // Cyan sphere
        Sphere(Vec3(1, 0, -2), 0.5, Vec3(1, 0, 1), 0.2),  // Magenta sphere
        Sphere(Vec3(0, -1, -2), 0.5, Vec3(1, 0.5, 0.5), 0.4) // Light red sphere
    };

    // Define the plane
    Plane plane(Vec3(0, -1, 0), Vec3(0, 1, 0)); // Plane at y = -1

    for (int frame = 0; frame < frame_count; ++frame) {
        auto start = std::chrono::high_resolution_clock::now();

        double angle = 2 * M_PI * frame / frame_count;
        Vec3 light_position(std::cos(angle) * radius, 1, std::sin(angle) * radius);

        std::ostringstream oss;
        oss << frame;
        std::string filename = "frame" + oss.str() + ".ppm";

        render_frame(filename, image_width, image_height, light_position, spheres, plane);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        std::cout << "Frame " << frame << " rendered in " << elapsed.count() << " seconds.\n";
    }

    return 0;
}

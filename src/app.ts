// https://raytracing.github.io/books/RayTracingInOneWeekend.html

const infinity = Number.POSITIVE_INFINITY;
const pi = Math.PI;

const degrees_to_radians = (degrees: number) => (degrees * pi) / 180;

// Returns a random real in [0,1).
const random_double = () => Math.random() % 1;

// Returns a random real in [min,max).
const random_double2 = (min: number, max: number) =>
    min + (max - min) * (Math.random() % 1);

const clamp = (x: number, min: number, max: number) =>
    Math.min(max, Math.max(min, x));

module Vec {
    export class Vec3 {
        constructor(public x: number, public y: number, public z: number) {}
        public neg() {
            return new Vec3(-this.x, -this.y, -this.z);
        }
        public add(vec: Vec3) {
            return new Vec3(this.x + vec.x, this.y + vec.y, this.z + vec.z);
        }
        public sub(vec: Vec3) {
            return new Vec3(this.x - vec.x, this.y - vec.y, this.z - vec.z);
        }
        public mutable_add(vec: Vec3) {
            this.x += vec.x;
            this.y += vec.y;
            this.z += vec.z;
        }
        public mul(vec: Vec3) {
            return new Vec3(this.x * vec.x, this.y * vec.y, this.z * vec.z);
        }
        public lengthSquared() {
            return this.x * this.x + this.y * this.y + this.z * this.z;
        }
        public length() {
            return Math.sqrt(this.lengthSquared());
        }
        public copy(vec: Vec3) {
            this.x = vec.x;
            this.y = vec.y;
            this.z = vec.z;
        }
        public near_zero() {
            // Return true if the vector is close to zero in all dimensions.
            return (
                Math.abs(this.x) < 1e-8 &&
                Math.abs(this.y) < 1e-8 &&
                Math.abs(this.z) < 1e-8
            );
        }
        public toString() {
            return "[" + this.x + " " + this.y + " " + this.z + "]";
        }
    }
}

import Vec3 = Vec.Vec3;
const vec3 = (x = 0, y = 0, z = 0) => new Vec3(x, y, z);
const zeroVec3 = vec3(0, 0, 0);

const add = (u: Vec3, v: Vec3) => vec3(u.x + v.x, u.y + v.y, u.z + v.z);
const sub = (u: Vec3, v: Vec3) => vec3(u.x - v.x, u.y - v.y, u.z - v.z);
const mul = (u: Vec3, v: Vec3) => vec3(u.x * v.x, u.y * v.y, u.z * v.z);
const cmul = (t: number, v: Vec3) => vec3(t * v.x, t * v.y, t * v.z);
const mulc = (v: Vec3, t: number) => vec3(t * v.x, t * v.y, t * v.z);
const divc = (v: Vec3, t: number) => vec3(v.x / t, v.y / t, v.z / t);

const unitVector = (v: Vec3) => divc(v, v.length());
const dot = (u: Vec3, v: Vec3) => u.x * v.x + u.y * v.y + u.z * v.z;
const cross = (u: Vec3, v: Vec3) =>
    vec3(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);

const lerp = (u: Vec3, v: Vec3, t: number) => cmul(1.0 - t, u).add(cmul(t, v));
const randomv = () => vec3(random_double(), random_double(), random_double());
const randomv2 = (min: number, max: number) =>
    vec3(
        random_double2(min, max),
        random_double2(min, max),
        random_double2(min, max)
    );
const reflect = (v: Vec3, n: Vec3) => v.sub(cmul(2 * dot(v, n), n));
const refract = (uv: Vec3, n: Vec3, etai_over_etat: number) => {
    const cos_theta = Math.min(dot(uv.neg(), n), 1);
    const r_out_perp = cmul(etai_over_etat, uv.add(cmul(cos_theta, n)));
    const r_out_parallel = cmul(
        -Math.sqrt(Math.abs(1.0 - r_out_perp.lengthSquared())),
        n
    );
    return r_out_perp.add(r_out_parallel);
};

import Point3 = Vec.Vec3;
const point3 = (x: number, y: number, z: number) => new Point3(x, y, z);
const zeroPoint3 = point3(0, 0, 0);

import Color = Vec.Vec3;
const color = (x: number, y: number, z: number) => new Color(x, y, z);
const zeroColor = color(0, 0, 0);

class Ray {
    constructor(public origin: Point3, public direction: Vec3) {}
    public at(t: number) {
        return this.origin.add(cmul(t, this.direction));
    }
    public copy(ray: Ray) {
        this.origin = ray.origin;
        this.direction = ray.direction;
    }
}
const ray = (orig: Point3, direction: Vec3) => new Ray(orig, direction);

const random_in_unit_disk = () => {
    while (true) {
        const p = vec3(random_double2(-1, 1), random_double2(-1, 1), 0);
        if (p.lengthSquared() < 1) return p;
    }
};

class Camera {
    public origin: Point3;
    public lower_left_corner: Point3;
    public horizontal: Vec3;
    public vertical: Vec3;
    public u: Vec3;
    public v: Vec3;
    public w: Vec3;
    public lens_radius: number;

    constructor(
        lookfrom: Point3,
        lookat: Point3,
        vup: Vec3,
        vfov: number, // vertical field-of-view in degrees
        aspect_ratio: number,
        aperture: number,
        focus_dist: number
    ) {
        const theta = degrees_to_radians(vfov);
        const h = Math.tan(theta / 2);
        const viewport_height = 2.0 * h;
        const viewport_width = aspect_ratio * viewport_height;

        this.w = unitVector(lookfrom.sub(lookat));
        this.u = unitVector(cross(vup, this.w));
        this.v = cross(this.w, this.u);

        this.origin = lookfrom;
        this.horizontal = cmul(focus_dist * viewport_width, this.u);
        this.vertical = cmul(focus_dist * viewport_height, this.v);
        this.lower_left_corner = this.origin
            .sub(divc(this.horizontal, 2))
            .sub(divc(this.vertical, 2))
            .sub(cmul(focus_dist, this.w));
        this.lens_radius = aperture / 2;
    }

    public get_ray(s: number, t: number) {
        const rd = cmul(this.lens_radius, random_in_unit_disk());
        const offset = mulc(this.u, rd.x).add(mulc(this.v, rd.y));
        return ray(
            this.origin.add(offset),
            this.lower_left_corner
                .add(cmul(s, this.horizontal))
                .add(cmul(t, this.vertical))
                .sub(this.origin)
                .sub(offset)
        );
    }
}

interface Material {
    scatter: (
        r: Ray,
        rec: HitRecord,
        attenuation: Color,
        scattered: Ray
    ) => boolean;
}

class HitRecord {
    public p!: Point3;
    public normal!: Vec3;
    public material?: Material;
    public t!: number;
    public front_face!: boolean;

    public set_face_normal(r: Ray, outward_normal: Vec3) {
        this.front_face = dot(r.direction, outward_normal) < 0;
        this.normal = this.front_face ? outward_normal : outward_normal.neg();
    }
    public copy(rec: HitRecord) {
        this.p = rec.p;
        this.normal = rec.normal;
        this.material = rec.material;
        this.t = rec.t;
        this.front_face = rec.front_face;
    }
}

interface Hittable {
    hit: (r: Ray, t_min: number, t_max: number, rec: HitRecord) => boolean;
}

class HittableList implements Hittable {
    public hittables: Hittable[] = [];

    constructor(hittable?: Hittable) {
        if (hittable) this.add(hittable);
    }

    public clear() {
        this.hittables.length = 0;
    }

    public add(hittable: Hittable) {
        this.hittables.push(hittable);
    }

    public hit(r: Ray, t_min: number, t_max: number, rec: HitRecord) {
        let temp_rec = new HitRecord();
        let hit_anything = false;
        let closest_so_far = t_max;

        for (const hittable of this.hittables) {
            if (hittable.hit(r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec.copy(temp_rec);
            }
        }
        return hit_anything;
    }
}

class Sphere implements Hittable {
    constructor(
        public center: Point3,
        public radius: number,
        public material: Material
    ) {}

    public hit(r: Ray, t_min: number, t_max: number, rec: HitRecord) {
        const oc = r.origin.sub(this.center);
        const a = r.direction.lengthSquared();
        const half_b = dot(oc, r.direction);
        const c = oc.lengthSquared() - this.radius * this.radius;
        const discriminant = half_b * half_b - a * c;

        if (discriminant < 0) {
            return false;
        }

        const sqrtd = Math.sqrt(discriminant);

        // Find the nearest root that lies in the acceptable range.
        let root = (-half_b - sqrtd) / a;
        if (root < t_min || t_max < root) {
            root = (-half_b + sqrtd) / a;
            if (root < t_min || t_max < root) return false;
        }

        rec.t = root;
        rec.p = r.at(root);
        const outward_normal = divc(rec.p.sub(this.center), this.radius);
        rec.set_face_normal(r, outward_normal);
        rec.material = this.material;

        return true;
    }
}

type DiffuseMethod = (normal?: Vec3) => Vec3;

const random_in_unit_sphere: DiffuseMethod = () => {
    while (true) {
        const p = randomv2(-1, 1);
        if (p.lengthSquared() < 1) return p;
    }
};

const random_unit_vector: DiffuseMethod = () =>
    unitVector(random_in_unit_sphere());

const random_in_hemisphere: DiffuseMethod = (normal?: Vec3) => {
    const in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal!) > 0.0)
        // In the same hemisphere as the normal
        return in_unit_sphere;
    else return in_unit_sphere.neg();
};

const diffuse_method: DiffuseMethod = random_unit_vector;

class Lambertian implements Material {
    constructor(public albedo: Color) {}

    public scatter(r: Ray, rec: HitRecord, attenuation: Color, scattered: Ray) {
        let scatter_direction = rec.normal.add(random_unit_vector());

        // Catch degenerate scatter direction
        if (scatter_direction.near_zero()) scatter_direction = rec.normal;

        scattered.copy(ray(rec.p, scatter_direction));
        attenuation.copy(this.albedo);
        return true;
    }
}

class Metal implements Material {
    public fuzz;

    constructor(public albedo: Color, fuzz: number) {
        this.fuzz = Math.min(fuzz, 1);
    }

    public scatter(
        r_in: Ray,
        rec: HitRecord,
        attenuation: Color,
        scattered: Ray
    ) {
        const reflected = reflect(unitVector(r_in.direction), rec.normal);
        scattered.copy(
            ray(rec.p, reflected.add(cmul(this.fuzz, random_in_unit_sphere())))
        );
        attenuation.copy(this.albedo);
        return dot(scattered.direction, rec.normal) > 0;
    }
}

class Dielectric implements Material {
    // ir = index_of_refraction
    constructor(public ir: number) {}

    private static reflectance(cosine: number, ref_idx: number) {
        // Use Schlick's approximation for reflectance.
        let r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * Math.pow(1 - cosine, 5);
    }

    public scatter(
        r_in: Ray,
        rec: HitRecord,
        attenuation: Color,
        scattered: Ray
    ) {
        attenuation.copy(color(1.0, 1.0, 1.0));
        const refraction_ratio = rec.front_face ? 1.0 / this.ir : this.ir;
        const unit_direction = unitVector(r_in.direction);

        const cos_theta = Math.min(dot(unit_direction.neg(), rec.normal), 1.0);
        const sin_theta = Math.sqrt(1.0 - cos_theta * cos_theta);

        const cannot_refract = refraction_ratio * sin_theta > 1.0;
        let direction: Vec3;

        if (
            cannot_refract ||
            Dielectric.reflectance(cos_theta, refraction_ratio) >
                random_double()
        )
            direction = reflect(unit_direction, rec.normal);
        else direction = refract(unit_direction, rec.normal, refraction_ratio);

        scattered.copy(ray(rec.p, direction));
        return true;
    }
}

const ray_color = (r: Ray, world: Hittable, depth: number): Color => {
    if (depth <= 0) return zeroColor;

    let rec = new HitRecord();
    if (world.hit(r, 0.00001, infinity, rec)) {
        let scattered = ray(zeroPoint3, zeroVec3);
        let attenuation = color(0, 0, 0);
        if (rec.material?.scatter(r, rec, attenuation, scattered)) {
            return attenuation.mul(ray_color(scattered, world, depth - 1));
        }
        return zeroColor;
    }
    const unit_direction = unitVector(r.direction);
    const t = 0.5 * (unit_direction.y + 1.0);
    return lerp(color(1.0, 1.0, 1.0), color(0.5, 0.7, 1.0), t);
};

const world1 = () => {
    const material_ground = new Lambertian(color(0.8, 0.8, 0.0));
    const material_center = new Lambertian(color(0.1, 0.2, 0.5));
    const material_left = new Dielectric(1.5);
    const material_right = new Metal(color(0.8, 0.6, 0.2), 0);

    const world = new HittableList();
    world.add(new Sphere(point3(0, -100.5, -1), 100, material_ground));
    world.add(new Sphere(point3(0, 0, -1), 0.5, material_center));
    world.add(new Sphere(point3(-1, 0, -1), 0.5, material_left));
    world.add(new Sphere(point3(-1, 0, -1), -0.45, material_left));
    world.add(new Sphere(point3(1, 0, -1), 0.5, material_right));

    return world;
};

const world2 = () => {
    const R = Math.cos(pi / 4);

    const material_left = new Lambertian(color(0, 0, 1));
    const material_right = new Lambertian(color(1, 0, 0));

    const world = new HittableList();
    world.add(new Sphere(point3(-R, 0, -1), R, material_left));
    world.add(new Sphere(point3(R, 0, -1), R, material_right));

    return world;
};

const random_scene = () => {
    const world = new HittableList();

    const ground_material = new Lambertian(color(0.5, 0.5, 0.5));
    world.add(new Sphere(point3(0, -1000, 0), 1000, ground_material));

    for (let a = -11; a < 11; a++) {
        for (let b = -11; b < 11; b++) {
            const choose_mat = random_double();
            const center = point3(
                a + 0.9 * random_double(),
                0.2,
                b + 0.9 * random_double()
            );

            if (center.sub(point3(4, 0.2, 0)).length() > 0.9) {
                let sphere_material: Material;

                if (choose_mat < 0.8) {
                    // diffuse
                    const albedo = randomv().mul(randomv());
                    const sphere_material = new Lambertian(albedo);
                    world.add(new Sphere(center, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    const albedo = randomv2(0.5, 1);
                    const fuzz = random_double2(0, 0.5);
                    const sphere_material = new Metal(albedo, fuzz);
                    world.add(new Sphere(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = new Dielectric(1.5);
                    world.add(new Sphere(center, 0.2, sphere_material));
                }
            }
        }
    }

    const material1 = new Dielectric(1.5);
    world.add(new Sphere(point3(0, 1, 0), 1.0, material1));

    const material2 = new Lambertian(color(0.4, 0.2, 0.1));
    world.add(new Sphere(point3(-4, 1, 0), 1.0, material2));

    const material3 = new Metal(color(0.7, 0.6, 0.5), 0.0);
    world.add(new Sphere(point3(4, 1, 0), 1.0, material3));

    return world;
};

const cam0 = (aspect_ratio: number) =>
    new Camera(
        point3(0, 0, 1),
        point3(0, 0, 0),
        vec3(0, 1, 0),
        90,
        aspect_ratio,
        0,
        1
    );

const cam1 = (aspect_ratio: number) =>
    new Camera(
        point3(-2, 2, 1),
        point3(0, 0, -1),
        vec3(0, 1, 0),
        90,
        aspect_ratio,
        0,
        1
    );

const cam2 = (aspect_ratio: number) =>
    new Camera(
        point3(-2, 2, 1),
        point3(0, 0, -1),
        vec3(0, 1, 0),
        20,
        aspect_ratio,
        0,
        1
    );

const cam3 = (aspect_ratio: number) => {
    const lookfrom = point3(3, 3, 2);
    const lookat = point3(0, 0, -1);
    const vup = vec3(0, 1, 0);
    const dist_to_focus = lookfrom.sub(lookat).length();
    const aperture = 2.0;
    return new Camera(
        lookfrom,
        lookat,
        vup,
        20,
        aspect_ratio,
        aperture,
        dist_to_focus
    );
};

const render = (
    ctx: CanvasRenderingContext2D,
    imageWidth: number,
    imageHeight: number,
    samples_per_pixel: number,
    max_depth: number,
    world: HittableList,
    cam: Camera
) => {
    const imageData = ctx.createImageData(imageWidth, imageHeight);

    let i = 0;
    for (let y = imageHeight; y >= 0; --y) {
        for (let x = 0; x < imageWidth; ++x) {
            const pixel_color = color(0, 0, 0);
            for (let s = 0; s < samples_per_pixel; ++s) {
                const u = (x + random_double()) / (imageWidth - 1);
                const v = (y + random_double()) / (imageHeight - 1);
                const r = cam.get_ray(u, v);
                pixel_color.mutable_add(ray_color(r, world, max_depth));
            }
            imageData.data[i] =
                Math.sqrt(pixel_color.x / samples_per_pixel) * 255;
            imageData.data[i + 1] =
                Math.sqrt(pixel_color.y / samples_per_pixel) * 255;
            imageData.data[i + 2] =
                Math.sqrt(pixel_color.z / samples_per_pixel) * 255;
            imageData.data[i + 3] = 255;
            i += 4;
        }
    }
    ctx.putImageData(imageData, 0, 0);
};

const main = () => {
    const run1 = (ctx: CanvasRenderingContext2D) => {
        const aspect_ratio = 16.0 / 9.0;
        const imageWidth = 300;
        const imageHeight = Math.floor(imageWidth / aspect_ratio);
        const samples_per_pixel = 20;
        const max_depth = 10;
        const world = world1();
        const cam = cam3(aspect_ratio);
        render(
            ctx,
            imageWidth,
            imageHeight,
            samples_per_pixel,
            max_depth,
            world,
            cam
        );
    };

    const run2 = (ctx: CanvasRenderingContext2D) => {
        const aspect_ratio = 3.0 / 2.0;
        const imageWidth = 400;
        const imageHeight = Math.floor(imageWidth / aspect_ratio);
        const samples_per_pixel = 5;
        const max_depth = 5;

        const world = random_scene();

        const lookfrom = point3(13, 2, 3);
        const lookat = point3(0, 0, 0);
        const vup = vec3(0, 1, 0);
        const dist_to_focus = 10;
        const aperture = 0.1;
        const cam = new Camera(
            lookfrom,
            lookat,
            vup,
            20,
            aspect_ratio,
            aperture,
            dist_to_focus
        );

        render(
            ctx,
            imageWidth,
            imageHeight,
            samples_per_pixel,
            max_depth,
            world,
            cam
        );
    };

    const canvas = document.getElementById(
        "canvas"
    ) as HTMLCanvasElement | null;
    const ctx = canvas?.getContext("2d");
    if (ctx) {
        run2(ctx);
    }
};

main();

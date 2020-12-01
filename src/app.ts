// https://raytracing.github.io/books/RayTracingInOneWeekend.html

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
        public mul(vec: Vec3) {
            return new Vec3(this.x * vec.x, this.y * vec.y, this.z * vec.z);
        }
        public lengthSquared() {
            return this.x * this.x + this.y * this.y + this.z * this.z;
        }
        public length() {
            return Math.sqrt(this.lengthSquared());
        }
        public toString() {
            return "[" + this.x + " " + this.y + " " + this.z + "]";
        }
    }
}

import Vec3 = Vec.Vec3;

const vec3 = (x = 0, y = 0, z = 0) => new Vec3(x, y, z);

const add = (u: Vec3, v: Vec3) => vec3(u.x + v.x, u.y + v.y, u.z + v.z);
const sub = (u: Vec3, v: Vec3) => vec3(u.x - v.x, u.y - v.y, u.z - v.z);
const mul = (u: Vec3, v: Vec3) => vec3(u.x * v.x, u.y * v.y, u.z * v.z);
const cmul = (t: number, v: Vec3) => vec3(t * v.x, t * v.y, t * v.z);
const mulc = (v: Vec3, t: number) => vec3(t * v.x, t * v.y, t * v.z);
const divc = (v: Vec3, t: number) => vec3(v.x / t, v.y / t, v.z / t);
const dot = (u: Vec3, v: Vec3) => u.x * v.x + u.y * v.y + u.z * v.z;
const cross = (u: Vec3, v: Vec3) =>
    vec3(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
const lerp = (u: Vec3, v: Vec3, t: number) => cmul(1.0 - t, u).add(cmul(t, v));

const unitVector = (v: Vec3) => divc(v, v.length());

import Point3 = Vec.Vec3;
const point3 = (x: number, y: number, z: number) => new Point3(x, y, z);

import Color = Vec.Vec3;
const color = (x: number, y: number, z: number) => new Color(x, y, z);

class Ray {
    constructor(public origin: Point3, public direction: Vec3) {}
    public at(t: number) {
        return this.origin.add(cmul(t, this.direction));
    }
}

const ray = (orig: Point3, direction: Vec3) => new Ray(orig, direction);

const hit_sphere = (center: Point3, radius: number, r: Ray) => {
    const oc = r.origin.sub(center);
    const a = dot(r.direction, r.direction);
    const b = 2 * dot(oc, r.direction);
    const c = dot(oc, oc) - radius * radius;
    const discriminant = b * b - 4 * a * c;
    return discriminant > 0;
};

const ray_color = (r: Ray): Color => {
    if (hit_sphere(point3(0, 0, -1), 0.5, r)) return color(1, 0, 0);

    const unit_direction = unitVector(r.direction);
    const t = 0.5 * (unit_direction.y + 1.0);
    return lerp(color(1.0, 1.0, 1.0), color(0.5, 0.7, 1.0), t);
};

const main = () => {
    const aspect_ratio = 16.0 / 9.0;
    const imageWidth = 400;
    const imageHeight = Math.floor(imageWidth / aspect_ratio);

    const viewport_height = 2.0;
    const viewport_width = aspect_ratio * viewport_height;
    const focal_length = 1.0;

    const origin = point3(0, 0, 0);
    const horizontal = vec3(viewport_width, 0, 0);
    const vertical = vec3(0, viewport_height, 0);
    const lower_left_corner = origin
        .sub(divc(horizontal, 2))
        .sub(divc(vertical, 2))
        .sub(vec3())
        .sub(vec3(0, 0, focal_length));

    const run = (ctx: CanvasRenderingContext2D) => {
        const imageData = ctx.createImageData(imageWidth, imageHeight);

        const render = () => {
            let i = 0;
            for (let y = imageHeight; y >= 0; --y) {
                for (let x = 0; x < imageWidth; ++x) {
                    const u = x / (imageWidth - 1);
                    const v = y / (imageHeight - 1);

                    const r = ray(
                        origin,
                        lower_left_corner
                            .add(cmul(u, horizontal))
                            .add(cmul(v, vertical))
                            .sub(origin)
                    );
                    const pixel_color = ray_color(r);

                    imageData.data[i] = pixel_color.x * 255;
                    imageData.data[i + 1] = pixel_color.y * 255;
                    imageData.data[i + 2] = pixel_color.z * 255;
                    imageData.data[i + 3] = 255;
                    i += 4;
                }
            }
            ctx.putImageData(imageData, 0, 0);
        };

        render();
    };

    const canvas = document.getElementById(
        "canvas"
    ) as HTMLCanvasElement | null;
    const ctx = canvas?.getContext("2d");
    if (ctx) {
        run(ctx);
    }
};

main();

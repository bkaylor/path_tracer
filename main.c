#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Converted from smallpt, a Path Tracer by Kevin Beason

#include "vec3.h"

// #define STDLIB_RAND

#ifndef STDLIB_RAND
#include "xoshiro256plus.c"
#include "splitmix64.c"
#endif

#define PI 3.14159265

typedef struct {
    vec3 origin, direction;
} Ray;

typedef enum {
    DIFFUSE, 
    SPECULAR, 
    REFRACTIVE
} Material;

typedef struct {
    float radius;
    vec3 position, emission, color;
    Material material;
} Sphere;

float random_float()
{
#ifdef STDLIB_RAND
    float random = rand() / (float)RAND_MAX;
    return random;
#else 
    return ((float)(xoshiro_next() >> 32)) / ((float)(0b11111111111111111111111111111111));
#endif
}

float sphere_intersect_with_ray(Sphere sphere, Ray ray)
{
    vec3 op = vec3_subtract(sphere.position, ray.origin);
    float t, epsilon = 1e-4;
    float b = vec3_dot_product(op, ray.direction);
    float determinant = b*b - vec3_dot_product(op, op) + sphere.radius*sphere.radius;

    if (determinant < 0) {
        return 0;
    } else {
        determinant = sqrt(determinant);
    }

    t = b - determinant;

    if (t > epsilon) {
        return t;
    } else {
        t = b + determinant;
        if (t > epsilon) {
            return t;
        } else {
            return 0;
        }
    }
}

// Global for scene data.
Sphere spheres[9];

float clamp(float f)
{
    if (f < 0) {
        return 0;
    } else if (f > 1) {
        return 1;
    } else {
        return f;
    }
}

int to_int(float f)
{
    return (int)(pow(clamp(f), 1/2.2) * 255 + 0.5);
}

bool intersect(Ray ray, float *t, int *id)
{
    int number_of_spheres = sizeof(spheres) / sizeof(Sphere);
    float d, infinity = *t = 1e20;

    for (int i = 0; i < number_of_spheres; i += 1)
    {
        d = sphere_intersect_with_ray(spheres[i], ray);
        if (d && d < *t) {
            *t = d;
            *id = i;
        }

    }

    return *t < infinity;
}

vec3 radiance(Ray ray, int depth)
{
    float distance_to_intersection = 0;
    int id = 0;

    if (!intersect(ray, &distance_to_intersection, &id)) return make_vec3(0, 0, 0);

    Sphere *hit_sphere = &spheres[id];

    if (depth > 10) return make_vec3(0, 0, 0);

    vec3 x = vec3_add(ray.origin, vec3_scalar_multiply(ray.direction, distance_to_intersection)); 
    vec3 n = vec3_normalize(vec3_subtract(x, hit_sphere->position));
    vec3 nl = vec3_dot_product(n, ray.direction) < 0 ? n : vec3_scalar_multiply(n, -1);
    vec3 f = hit_sphere->color;

    float p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;

    depth += 1;

    if (depth > 5) 
    {
        if (random_float() < p) 
        { 
            f = vec3_scalar_multiply(f, 1/p);
        } 
        else 
        {
            return hit_sphere->emission;
        }
    } 

    switch (hit_sphere->material)
    {
        case DIFFUSE: {
            /*
            double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2); 
            Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u; 
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 
            return obj.e + f.mult(radiance(Ray(x,d),depth,Xi)); 
            */
            float r1 = 2 * PI * random_float(); 
            float r2 = random_float(); 
            float r2s = sqrt(r2);

            vec3 w = nl; 
            vec3 u = vec3_normalize(vec3_cross_product(fabs(w.x) > 0.1f ? make_vec3(0, 1, 0) : make_vec3(1, 0, 0), w));
            vec3 v = vec3_cross_product(w, u);
            vec3 d = vec3_normalize(vec3_add3(vec3_scalar_multiply(u, cos(r1) * r2s),
                                              vec3_scalar_multiply(v, sin(r1) * r2s),
                                              vec3_scalar_multiply(w, sqrt(1 - r2))));
            Ray new_ray;
            // NOTE(bkaylor): Big thanks to @UglySwedishFisk on Discord for the next line suggestion. 
            new_ray.origin = vec3_add(x, vec3_scalar_multiply(nl, 0.03f));
            // new_ray.origin = x;
            new_ray.direction = d;
            return vec3_add(hit_sphere->emission, vec3_linear_multiply(radiance(new_ray, depth), f));
        } break;
        case SPECULAR: {
            /*
            return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi)); 
            */
            Ray new_ray;
            new_ray.origin = x;
            new_ray.direction = vec3_subtract(ray.direction, vec3_scalar_multiply(n, 2 * vec3_dot_product(n, ray.direction)));
            return vec3_add(hit_sphere->emission, vec3_linear_multiply(radiance(new_ray, depth), f));
        } break;
        case REFRACTIVE: {
           /*
           Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION 
           bool into = n.dot(nl)>0;                // Ray from outside going in? 
           double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t; 
           if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection 
             return obj.e + f.mult(radiance(reflRay,depth,Xi)); 
           Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
           double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
           double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P); 
           return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette 
             radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) : 
             radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr); 
           */
            Ray reflection_ray;
            reflection_ray.origin = x;
            reflection_ray.direction = vec3_subtract(ray.direction, vec3_scalar_multiply(n, 2 * vec3_dot_product(n, ray.direction)));

            bool into = vec3_dot_product(n, nl) > 0;
            float nc = 1;
            float nt = 1.5;
            float nnt = into ? nc / nt : nt / nc;
            float ddn = vec3_dot_product(ray.direction, nl);

            float cos2t = 1 - nnt * nnt * (1 - ddn * ddn);

            if (cos2t < 0)
            {
                return vec3_add(hit_sphere->emission, vec3_linear_multiply(radiance(reflection_ray, depth), f));
            }

            vec3 tdir = vec3_normalize(vec3_subtract(vec3_scalar_multiply(ray.direction, nnt),
                                                     vec3_scalar_multiply(n, (into ? 1 : -1) * ddn * nnt + sqrt(cos2t))));  
            float a = nt - nc;
            float b = nt + nc;
            float r0 = a*a / (b*b);
            float c = 1 - (into ? -ddn : vec3_dot_product(tdir, n));

            float Re = r0 + (1 - r0) * c * c * c * c * c;
            float Tr = 1 - Re;
            float P = 0.25 + 0.5 * Re;
            float RP = Re / P;
            float TP = Tr / (1 - P);

            return vec3_add(hit_sphere->emission, depth > 2 ? 
                                                      (random_float() < P ? 
                                                          vec3_scalar_multiply(radiance(reflection_ray, depth), RP) :
                                                          vec3_scalar_multiply(radiance((Ray){x, tdir}, depth), TP)) :
                                                      vec3_add(
                                                          vec3_scalar_multiply(radiance(reflection_ray, depth), Re),
                                                          vec3_scalar_multiply(radiance((Ray){x, tdir}, depth), Tr)));


        } break;
        default: {
            return make_vec3(0, 0, 0);
        }
    }

    return make_vec3(0, 0, 0);
}

int main(int argc, char *argv[])
{
    spheres[0] = (Sphere){1e5,  (vec3){ 1e5+1,40.8,81.6 }, (vec3){0,0,0      }, (vec3){.75,.25,.25   }, DIFFUSE};   //Left 
    spheres[1] = (Sphere){1e5,  (vec3){-1e5+99,40.8,81.6}, (vec3){0,0,0      }, (vec3){.25,.25,.75   }, DIFFUSE};   //Rght 
    spheres[2] = (Sphere){1e5,  (vec3){50,40.8, 1e5     }, (vec3){0,0,0      }, (vec3){.75,.75,.75   }, DIFFUSE};   //Back 
    spheres[3] = (Sphere){1e5,  (vec3){50,40.8,-1e5+170 }, (vec3){0,0,0      }, (vec3){0, 0, 0       }, DIFFUSE};   //Frnt 
    spheres[4] = (Sphere){1e5,  (vec3){50, 1e5, 81.6    }, (vec3){0,0,0      }, (vec3){.75,.75,.75   }, DIFFUSE};   //Botm 
    spheres[5] = (Sphere){1e5,  (vec3){50,-1e5+81.6,81.6}, (vec3){0,0,0      }, (vec3){.75,.75,.75   }, DIFFUSE};   //Top 
    spheres[6] = (Sphere){16.5, (vec3){27,16.5,47       }, (vec3){0,0,0      }, (vec3){.999,.999,.999}, SPECULAR};  //Mirr 
    spheres[7] = (Sphere){16.5, (vec3){73,16.5,78       }, (vec3){0,0,0      }, (vec3){.999,.999,.999}, REFRACTIVE};//Glas 
    spheres[8] = (Sphere){600,  (vec3){50,681.6-.27,81.6}, (vec3){6,6,6      }, (vec3){0, 0, 0       }, DIFFUSE};   //Lite 

#ifdef STDLIB_RAND
    srand(time(NULL));
#else
    splitmix_x = 8189475927166758637;
    s[0] = splitmix_next();
    s[1] = splitmix_next();
    s[2] = splitmix_next();
    s[3] = splitmix_next();
#endif

    int width = 2560, height = 1440;
    int samples = 1;
    if (argc == 2) samples = atoi(argv[1]) / 4;

    Ray camera;
    camera.origin = make_vec3(50, 52, 295.6);
    camera.direction = vec3_normalize(make_vec3(0, -0.042612, -1));

    vec3 x_increment = make_vec3((float)width * 0.5135 / (float)height, 0, 0);
    vec3 y_increment = vec3_scalar_multiply(vec3_normalize(vec3_cross_product(x_increment, camera.direction)), 0.5135f);

    vec3 r;

    vec3 *image = malloc(sizeof(vec3) * width * height);
    for(int i = 0; i < width * height; i += 1) image[i] = make_vec3(0, 0, 0);

    for (int y = 0; y < height; y += 1) // Rows
    {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", 
                samples * 4, 
                100.0f * y/(height - 1)); 

        for (int x = 0; x < width; x += 1) // Columns
        {
            for (int subpixel_y = 0, i = (height-y-1) * width + x; subpixel_y < 2; subpixel_y += 1) // Subpixel rows
            {
                for (int subpixel_x = 0; subpixel_x < 2; subpixel_x += 1) // Subpixel columns
                {
                    r = make_vec3(0, 0, 0);

                    for (int s = 0; s < samples; s += 1)
                    {
                        float r1 = 2*random_float();
                        float dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);

                        float r2 = 2*random_float();
                        float dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

                        vec3 d_first_term =  vec3_scalar_multiply(x_increment, (((subpixel_x + 0.5 + dx)/2 + x)/width - 0.5));
                        vec3 d_second_term = vec3_scalar_multiply(y_increment, (((subpixel_y + 0.5 + dy)/2 + y)/height - 0.5));

                        vec3 d = vec3_add3(d_first_term, d_second_term, camera.direction);

                        Ray ray;
                        ray.origin = vec3_add(camera.origin, vec3_scalar_multiply(d, 140.0f));
                        ray.direction = vec3_normalize(d);
                        r = vec3_add(r, vec3_scalar_multiply(radiance(ray, 0), 1.0f / samples));
                    }

                    image[i] = vec3_add(image[i], vec3_scalar_multiply(make_vec3(clamp(r.x), clamp(r.y), clamp(r.z)), 0.25f));
                }
            }
        }
    }

    FILE *f = fopen("image.ppm", "w");
    fprintf(f, "P3\n%d %d\n%d\n", width, height, 255); 
    for (int i = 0; i < width * height; i += 1) 
    {
        fprintf(f,"%d %d %d ", to_int(image[i].x), to_int(image[i].y), to_int(image[i].z)); 
    }
}

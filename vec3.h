typedef struct {
    float x;
    float y;
    float z;
} vec3;

static inline vec3 make_vec3(float x, float y, float z)
{
    return (vec3){x, y, z};
}

static inline vec3 vec3_add(vec3 a, vec3 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

static inline vec3 vec3_add3(vec3 a, vec3 b, vec3 c)
{
    return vec3_add(vec3_add(a, b), c);
}

static inline vec3 vec3_subtract(vec3 a, vec3 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

static inline vec3 vec3_scalar_multiply(vec3 a, float b)
{
    a.x *= b;
    a.y *= b;
    a.z *= b;
    return a;
}

static inline vec3 vec3_linear_multiply(vec3 a, vec3 b)
{
    vec3 temp;
    temp.x = a.x * b.x;
    temp.y = a.y * b.y;
    temp.z = a.z * b.z;
    return temp;
}

static inline vec3 vec3_cross_product(vec3 a, vec3 b)
{
    vec3 temp;
    temp.x = a.y * b.z - a.z * b.y;  
    temp.y = a.z * b.x - a.x * b.z;
    temp.z = a.x * b.y - a.y * b.x;
    return temp;
}

static inline float vec3_dot_product(vec3 a, vec3 b)
{
    vec3 temp;
    temp.x = a.x * b.x;
    temp.y = a.y * b.y;
    temp.z = a.z * b.z;
    return temp.x + temp.y + temp.z;
}

static inline vec3 vec3_normalize(vec3 a)
{
    float scale;
    scale = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);

    a.x /= scale;
    a.y /= scale;
    a.z /= scale;

    return a;
}

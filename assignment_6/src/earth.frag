//=============================================================================
//
//   Exercise code for the lecture "Introduction to Computer Graphics"
//     by Prof. Mario Botsch, Bielefeld University
//
//   Copyright (C) by Computer Graphics Group, Bielefeld University
//
//=============================================================================

#version 140

in vec3 v2f_normal;
in vec2 v2f_texcoord;
in vec3 v2f_light;
in vec3 v2f_view;

out vec4 f_color;

uniform sampler2D day_texture;
uniform sampler2D night_texture;
uniform sampler2D cloud_texture;
uniform sampler2D gloss_texture;
uniform bool greyscale;

const float shininess = 20.0;
const vec3  sunlight = vec3(1.0, 0.941, 0.898);

void main()
{
    /** \todo
    * - Copy your working code from the fragment shader of your Phong shader use it as
    * starting point
    * - instead of using a single texture, use the four texures `day_texure`, `night_texure`,
    * `cloud_texure` and `gloss_texture` and mix them for enhanced effects
    * Hints:
    * - cloud and gloss textures are just greyscales. So you'll just need one color-
    * component.
    * - The texture(texture, 2d_position) returns a 4-vector (rgba). You can use
    * `texture(...).r` to get just the red component or `texture(...).rgb` to get a vec3 color
    * value
    * - use mix(vec3 a,vec3 b, s) = a*(1-s) + b*s for linear interpolation of two colors
     */

    vec3 r = normalize(reflect(-v2f_light,v2f_normal));

    vec3 Ia = 0.2 * sunlight;
    float d_test = dot(v2f_normal, v2f_light);
    float d = max(0.0, d_test);
    float s = 0.0;
    if(d_test >= 0){
        s = max(0.0,pow(dot(r, v2f_view), shininess));
    }

    vec3 day = texture(day_texture, v2f_texcoord).rgb;
    float cloudiness = texture(cloud_texture, v2f_texcoord).r;
    vec3 night = texture(night_texture, v2f_texcoord).rgb * (1.0 - cloudiness);
    float gloss = texture(gloss_texture, v2f_texcoord).r * (1.0 - cloudiness);
    
    day = mix(Ia * day + d * sunlight * day + s * gloss * sunlight, Ia * vec3(cloudiness) + d * sunlight * vec3(cloudiness), cloudiness);
    
    vec3 color = mix(night, day, d);

    if (greyscale) color = vec3(0.299*color.r+0.587*color.g+0.114*color.b);

    f_color = vec4(color, 1.0);

}

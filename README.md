srpbr is a simple pure c++ 11 software-renderer with pbr shading model, easy to understand, ease to debug.

#feature
1. left-hand z-up coordinate system
2. software-rasterizer
3. perspective-correct interpolation
4. inverse z-buffer (z-test / z-write)
5. back culling / frustum culling / viewport scissor 
6. support wavefront .obj file
7. physically-based-shading with ibl

ibl precompute
generate ibl textures from an environment cubemap
1. indirect diffuse lighting
precompute indirect diffuse lighting stored in an irradiance map from a radiance map(environment cubemap) 
2. indirect specular lighting
generate_irradiance_map("./resource/ibl_textures/env.png", "./resource/ibl_textures/irradiance.png");
generate_prefilter_envmap("./resource/ibl_textures/env.png", "./resource/ibl_textures/prefilter");
generate_BRDF_LUT("./resource/brdf_lut.png");


|direct lighting only|
| ------------- |
|![](https://github.com/niepp/srpbr/blob/main/images/direct%20lighting%20only.png)|


|with ibl|
| ------------- |
|![](https://github.com/niepp/srpbr/blob/main/images/constant%20color%20with%20ibl.png)|


|no indirect lighting|
| ------------- |
|![](https://github.com/niepp/srpbr/blob/main/images/no%20indirect%20lighting.png)|

|no specular|
| ------------- |
|![](https://github.com/niepp/srpbr/blob/main/images/no%20specular.png)|

|standard pbr with ibl|
| ------------- |
|![](https://github.com/niepp/srpbr/blob/main/images/standard%20pbr%20with%20ibl.png)|



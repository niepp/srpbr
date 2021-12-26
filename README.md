# **srpbr** is a simple software-renderer with pbr shading model, easy to understand, easy to debug.

## features</br>

1. pure c++ 11 implementation
2. left-hand z-up coordinate system
3. software-rasterizer
4. perspective-correct interpolation
5. z-buffer (reverse-z tech / z-test / z-write)
6. back culling / frustum culling / viewport clip
7. support wavefront .obj file
8. physically-based-shading with ibl
9. bilinear / trilinear texture sample
10. cube map sample
11. reinhard toon-mapping


## ibl precompute
generate ibl textures from an environment cubemap
1. indirect diffuse lighting
    - precompute indirect diffuse lighting stored in an irradiance map from a radiance map(environment cubemap) 
    ```cpp
    generate_irradiance_map("./resource/ibl_textures/env.png", "./resource/ibl_textures/irradiance.png");
    ```
2. indirect specular lighting
    - Split Sum Approximation, result in two part
    ```cpp
    generate_prefilter_envmap("./resource/ibl_textures/env.png", "./resource/ibl_textures/prefilter");
    generate_BRDF_LUT("./resource/brdf_lut.png");
    ```



## result

|depth buffer|
| ------------- |
|![](https://github.com/niepp/srpbr/blob/main/result/depth_gun.png)|

|direct lighting only|
| ------------- |
|![](https://github.com/niepp/srpbr/blob/main/result/direct_lighting_only.jpg)|

|with ibl|
| ------------- |
|![](https://github.com/niepp/srpbr/blob/main/result/constant_color_with_ibl.jpg)|

|standard pbr with ibl|
| ------------- |
|![](https://github.com/niepp/srpbr/blob/main/result/srpbr.jpg)|

## see more

 https://niepp.github.io/2016/05/20/PBR-notes.html
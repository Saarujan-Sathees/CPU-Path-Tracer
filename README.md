# CPU-Path-Tracer
A simple static, but slow path tracer that uses the CPU to render 3D scenes. Contains the following:

- The option to render asynchronously, and the option to render to a PPM file or to the offline webpage image-viewer.html that allows frames to be rendered per aliasing count (rays-per-pixel)
- Currently supports spheres, quadrilaterals and triangles with the use of bounding-box volumes
- Supports field-of-view, focus distance & blur, different lens radii, and a moveable camera
- Supports limited material settings:
    - Colors
    - Specular / Diffuse Materials
    - Lights
- A Terrain object that is created using a variation of the Marching Cubes algorithm
- A custom-made JSON parser that allows scene settings and information to be stored as a JSON file

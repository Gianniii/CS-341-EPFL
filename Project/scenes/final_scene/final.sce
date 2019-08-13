# camera: eye-point, look-at-point, up, fovy, width, height
camera  -1.8 1.6 0.35  0 0 0.5  0 1 0  90  1920 1080 0.37

# recursion depth
depth  4

# background color
background 1 1 1

# global ambient light
ambience   0.1 0.1 0.1

# light: position and color
light -10 20 10  0.01 0.01 0.01
light -15 10 15  0.5 0.5 0.5
#light -6 6 6  0.1 0.1 0.1
# meshes: filename, shading, material (ambient, diffuse, specular, shininess)

mesh joint.obj  PHONG  0.0 0.0 0.0 0.0 0.0 0.0  0.0 0.0 0.0  3.0   0.0   0.0  0.0 0.0 0.0
mesh glass.obj PHONG  0.39 0.93 0.64  0.39 0.93 0.64   1.0 1.0 1.0  900.0   0.3   0.0  0.0 0.0 0.0
mesh wall1.obj PHONG  0.5 0.5 0.5 0.0 0.0 0.0  0.0 0.0 0.0  0.0   0.0   0.0  0.0 0.0 0.0
mesh sand.obj PHONG  0.95 0.87 0.60  0.95 0.87 0.60  1.0 1.0 1.0  15.0  0.0  0.0  0.0 0.0 0.0
mesh wall2.obj PHONG  0.5 0.5 0.5 0.0 0.0 0.0  0.0 0.0 0.0  3.0   0.0   0.0  0.0 0.0 0.0

sphere  -1.4 1.3 0.8  0.25  0.8 0.0 0.0  0.8 0.0 0.0  1.0 1.0 1.0  100.0  0.0  1.5  0.0 0.18 0.08

#sol
#plane  2 0 0  0 1 0  0.95 0.87 0.60  0.95 0.87 0.60  1.0 1.0 1.0  10.0  0.1  0.0  0.0 0.0 0.0

//=============================================================================
//
//   Exercise code for the lecture
//   "Introduction to Computer Graphics"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) Computer Graphics Group, Bielefeld University.
//
//=============================================================================

//== includes =================================================================

#include "StopWatch.h"
#include "Scene.h"

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

/// Program entry point.
int main(int argc, char **argv) {
    // Parse input scene file/output path from command line arguments
    struct RaytraceJob { std::string scenePath, outPath; };
    std::vector<RaytraceJob> jobs;

    if (argc == 3)
        jobs.emplace_back(RaytraceJob{argv[1], argv[2]});
    else if ((argc == 2) && argv[1][0] == '0') {
        jobs = { {
            {"../scenes/spheres/spheres.sce",       "spheres.png"},
            {"../scenes/cylinders/cylinders.sce",   "cylinders.png"},
            {"../scenes/combo/combo.sce",           "combo.png"},
            {"../scenes/molecule/molecule.sce",     "molecule.png"},
            {"../scenes/molecule2/molecule2.sce",   "molecule2.png"},
            {"../scenes/cube/cube.sce",             "cube.png"},
            {"../scenes/mask/mask.sce",             "mask.png"},
            {"../scenes/mirror/mirror.sce",         "mirror.png"},
            {"../scenes/toon_faces/toon_faces.sce", "toon_faces.png"},
            {"../scenes/office/office.sce",         "office.png"},
            {"../scenes/rings/rings.sce",           "rings.png"}
        } };
    }
    else {
        std::cerr << "Usage: " << argv[0] << " input.sce output.png\n";
        std::cerr << "Or: " << argv[0] << " 0\n";
        std::cerr << std::flush;
        exit(1);
    }

    for (const auto &job : jobs) {
        std::cout << "Read scene '" << job.scenePath << "'..." << std::flush;
        Scene s(job.scenePath);
        std::cout << "\ndone (" << s.numObjects() << " objects)\n";

        StopWatch timer;
        std::cout << "Ray tracing..." << std::flush;
        timer.start();
        auto image = s.render();
        timer.stop();
        std::cout << " done (" << timer << ")\n";

        std::cout << "Write image...";
        image.write(job.outPath);
        std::cout << "done\n";
    }
}

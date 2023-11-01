#include "Headers/Vector.h"
#include "Headers/Ray.h"
#include "Headers/Objects.h"
#include "Headers/Camera.h
#include "Headers/Frame.h"
#include "Headers/Scene.h"
#include "Headers/Random.h"
#include "Headers/Render.h"
#include "Headers/Parser.h"
#include <time.h>
#include <experimental/source_location>

void currTime() {
    std::time_t current = std::time(0);
    std::tm* time = std::localtime(&current);

    std::printf("%02d:%02d:%02d\n", time -> tm_hour > 12 ? time -> tm_hour - 12 : time -> tm_hour, time -> tm_min, time -> tm_sec);
}

void lapTime(std::clock_t& start, std::experimental::source_location curr = std::experimental::source_location::current()) {
    std::clock_t end = std::clock();
    std::cout << "Line " << curr.line() << ": " << double(end - start) / CLOCKS_PER_SEC << "s\n";
    start = end; 
}

int main() {
    Camera* camera = new Camera();
    Scene* scene = new Scene();
    std::clock_t begin;
    currTime();

    if (std::filesystem::exists("Settings.json")) {
        begin = std::clock();
        JSON settings = Parser::parseJSON("Settings.json");
        Types* renderType = settings.get("frame_update_on");

        if (renderType == nullptr || renderType -> boolValue == false) { //PPM Image
            Frame* frame = Parser::buildFrame(*camera, *scene, settings);
            frame -> render(camera, scene);
        } else { //Frame Viewer
            FrameServer* frame = Parser::buildFrameServer(*camera, *scene, settings);
            frame -> render(camera, scene);
        }

    } else {
        short renderType;
        std::cout << "1 - Single Image\n2 - Frame Update\n\nChoose: ";
        std::cin >> renderType;
        std::cout << "\n";
        begin = std::clock();
        
        scene -> push(new Sphere(Vector(20, 160, -1), 80, new Light(Color(255, 255, 255), 10)));
        //scene -> push(new Sphere(Vector(0, 0, -1), 0.5, new Material(Color(50, 100, 250))));
        scene -> push(new Triangle(Vector(5, 1, -1), Vector(0, 0, -1), Vector(-2, 2, -1), new Material(Color(255, 0, 0))));
        scene -> push(new Sphere(Vector(0, -100.5, -1), 100, new Material(Color(100, 100, 100))));

        if (renderType == 1) {
            Frame* frame = new Frame(1280, 720, 100, "output.ppm");
            camera = new Camera(Vector(-2, 2, 1), Vector(0, 0, -1), 90, frame -> getAspectRatio());
            frame -> render(camera, scene);
        } else {
            FrameServer* frame = new FrameServer(1280, 720, 1, 2021);
            Camera* camera = new Camera(Vector(-2, 2, 1), Vector(0, 0, -1), 90, frame -> getAspectRatio());
            frame -> render(camera, scene);
        }
    }
    
    lapTime(begin);
    currTime();
}
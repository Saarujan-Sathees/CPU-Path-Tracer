#ifndef PARSER
#define PARSER

#include "Vector.h"
#include "Frame.h"
#include "Objects.h"
#include "Camera.h"
#include "Materials.h"
#include "Terrain Generation/Terrain.h"
#include "JSON.h"
#include <string>
#include <fstream>

class Parser {
    private:
    static bool isDigit(char c) {
        return c >= '0' && c <= '9';
    }

    static Types* parseBoolean(std::string& line) {
        Types* res = new Types();
        res -> boolValue = line[0] == 't';
        
        line = line.substr(5 - res -> boolValue); //4 for true, 5 for false
        return res;
    }

    static Types* parseNumber(std::string& line) {
        int index = 0;
        bool isDecimal = false, isPositive = line[index] != '-';
        Types* value = new Types();
        value -> dblValue = isPositive ? line[index] - '0' : 0;

        double decMult = 0.1;
    
        for (index = index + 1; index < line.size(); ++index) {
            if (isDigit(line[index])) { 
                if (isDecimal) {
                    value -> dblValue += (line[index] - '0') * decMult;
                    decMult *= 0.1;
                } else {
                    value -> dblValue = value -> dblValue * 10 + (line[index] - '0');
                }
            } else if (!isDecimal && line[index] == '.') {
                isDecimal = true;
            } else {
                value -> dblValue *= isPositive * 2 - 1; //Converting to negative if a '-' was read
                line = line.substr(index);
                return value;
            }
        }

        value -> dblValue *= isPositive * 2 - 1; //Converting to negative if a '-' was read
        line = line.substr(index);
        return value;
    } 

    static void parseSpace(std::string& line) {
        int index = 0;
        while (line[index] == ' ' || line[index] == '\t' || line[index] == '\n') {
            ++index;
        }

        line = line.substr(index);
    }

    static void parseComment(std::string& file) {
        int index = 0;
        while (file[index] != '/' || file[index - 1] != '*') {
            ++index;
        }

        file = file.substr(index + 1);
    }

    static std::string parseKey(std::string& line) {
        std::string key = "";
        int index = 0;
        while (line[index] != '"') {
            ++index;
        }
        
        while (line[++index] != '"') {
            key += (char) tolower((int) line[index]);
        }

        line = line.substr(index + 1);
        return key;
    }

    static Types* parseVector(std::string& line) {
        std::string strValue = parseKey(line);
        Types* value = new Types();
        double xyz[3];

        for (int i = 0; i < 3; ++i) {
            xyz[i] = parseNumber(strValue) -> dblValue;
            parseSpace(strValue);
        }

        value -> vecValue = new Vector(xyz[0], xyz[1], xyz[2]);
        return value;
    }

    static Types* parseArray(std::string& file) {
        Types* value = new Types();
        file = file.substr(1); //Removing the open square bracket
        while (file[0] != ']') {
            if (file[0] == '/' && file[1] == '*') { //Comment
                parseComment(file);
            } else if (file[0] == ',') { //Comma
                file = file.substr(1);
            } else if (file[0] == ' ' || file[0] == '\n' || file[0] == '\t') { //Space
                parseSpace(file);
            } else if (file[0] == 't' || file[0] == 'f') {
                value -> arrValue.push_back(parseBoolean(file));
            } else if (file[0] == '-' || isDigit(file[0]) || file[0] == '.') { //Number
                value -> arrValue.push_back(parseNumber(file));
            } else if (file[0] == '{') { //Nested JSON
                value -> arrValue.push_back(new Types(parse(file)));
            } else if (file[0] == '"') { //Vector or File Path
                if (file[1] == '.' && file[2] == '/')  //File Path
                    value -> arrValue.push_back(new Types(parseKey(file)));
                else //Vector
                    value -> arrValue.push_back(parseVector(file));
            } else if (file[0] == '[') { //Array
                value -> arrValue.push_back(parseArray(file));
            } 
        }

        file = file.substr(1);
        return value;
    }

    static JSON parse(std::string& file) {
        int index = 0;
        while (file[index] != '{') {
            ++index;
        }

        file = file.substr(index + 1);
        JSON res = JSON();
        std::string currKey;
        Types* currValue;

        while (file[0] != '}') {
            if (file[0] == '/' && file[1] == '*') { //Comment
                parseComment(file);
            } else if (file[0] == ':' || file[0] == ',') { //Other Characters
                file = file.substr(1);
            } else if (file[0] == ' ' || file[0] == '\n' || file[0] == '\t') { //Space
                parseSpace(file);
            } else if (file[0] == 't' || file[0] == 'f') { //Boolean
                res.addPair(currKey, parseBoolean(file));
                currKey = "";
            } else if (file[0] == '-' || isDigit(file[0]) || file[0] == '.') { //Number
                res.addPair(currKey, parseNumber(file));
                currKey = "";
            } else if (file[0] == '{') { //Nested JSON
                res.addPair(currKey, new Types(parse(file)));
                currKey = "";
            } else if (file[0] == '"') { //Vector, File Path or Key
                if (currKey == "") { //Key
                    currKey = parseKey(file);
                } else if (file[1] == '.' && file[2] == '/') { //File Path
                    res.addPair(currKey, new Types(parseKey(file)));
                    currKey = "";
                } else { //Vector
                    res.addPair(currKey, parseVector(file));
                    currKey = "";
                }
            } else if (file[0] == '[') { //Array
                res.addPair(currKey, parseArray(file));
                currKey = "";
            } 
        }

        file = file.substr(1);
        return res;
    }

    public:
    static JSON parseJSON(const std::string& fileName) {
        std::ifstream input(fileName);
        std::string file = "";
        
        while (!input.eof()) {
            file += (char) input.get();
        }
        
        return parse(file);
    }

    static Scene* buildScene(const JSON& json, const Camera& camera) {
        JSON objectList = json.get("objects") -> jsonValue;
        Scene* scene = new Scene();

        if (objectList.contains("spheres")) {
            for (Types* obj : objectList.get("spheres") -> arrValue) {
                Traceable* toAdd = new Sphere();
                toAdd -> buildObject(obj -> jsonValue);
                scene -> push(toAdd);
            }
        }

        if (objectList.contains("triangles")) {
            for (Types* obj : objectList.get("triangles") -> arrValue) {
                Triangle* toAdd = new Triangle();
                toAdd -> buildObject(obj -> jsonValue);
                toAdd -> correctNormal(camera.origin);
                scene -> push(toAdd);
            }
        }

        if (objectList.contains("rectangles")) {
            for (Types* obj : objectList.get("rectangles") -> arrValue) {
                Quadrilateral* toAdd = new Quadrilateral();
                toAdd -> buildObject(obj -> jsonValue);
                toAdd -> correctNormal(camera.origin);
                scene -> push(toAdd);
            }
        }

        if (objectList.contains("aligned_rectangles")) {
            for (Types* obj : objectList.get("aligned_rectangles") -> arrValue) {
                Traceable* toAdd = new AlignedRectangle();
                toAdd -> buildObject(obj -> jsonValue);
                scene -> push(toAdd);
            }
        }

        if (objectList.contains("terrains")) {
            for (Types* obj : objectList.get("terrains") -> arrValue) {
                Terrain* toAdd = new Terrain();
                toAdd -> buildObject(obj -> jsonValue);
                toAdd -> Generate(camera.origin);
                scene -> push(toAdd);
            }
        }

        return scene;
    }

    static Camera* buildCamera(const JSON& json) {
        Vector position = Vector(0, 0, 0), target = Vector(0, 0, -1);
        short fieldOfView = 90;
        double focusDist = -1, aperture = 0, aspectRatio;

        Types* temp = json.get("camera_position");
        if (temp != nullptr) 
            position = *(temp -> vecValue);
        
        temp = json.get("camera_target");
        if (temp != nullptr) 
            target = *(temp -> vecValue);
        
        temp = json.get("field_of_view");
        if (temp != nullptr) 
            fieldOfView = (short) temp -> dblValue;
        
        temp = json.get("aperture");
        if (temp != nullptr) 
            aperture = temp -> dblValue;
            
        temp = json.get("focus_distance");
        if (temp != nullptr) 
            focusDist = temp -> dblValue;
        
        temp = json.get("frame_width");
        if (temp != nullptr) 
            aspectRatio = temp -> dblValue;

        temp = json.get("frame_height");
        if (temp != nullptr) 
            aspectRatio /= temp -> dblValue;

        return new Camera(position, target, fieldOfView, aspectRatio, focusDist, aperture);
    }

    static Frame* buildFrame(Camera& camera, Scene& scene, const JSON& json = parseJSON("Settings.json")) {
        int width, height;
        short raysPerPixel, threadCount = 1;
        Vector backgroundColor = BLACK;

        std::string filePath = "output.ppm";
        Types* temp = json.get("frame_width");
        if (temp != nullptr) 
            width = (int) temp -> dblValue;

        temp = json.get("frame_height");
        if (temp != nullptr) 
            height = (int) temp -> dblValue;

        temp = json.get("rays_per_pixel");
        if (temp != nullptr) 
            raysPerPixel = (short) temp -> dblValue;

        temp = json.get("output_file");
        if (temp != nullptr) 
            filePath = temp -> strValue;

        temp = json.get("background_color");
        if (temp != nullptr)
            backgroundColor = Color(*(temp -> vecValue));

        temp = json.get("thread_count");
        if (temp != nullptr)
            threadCount = (short) temp -> dblValue;

        camera = *buildCamera(json);
        scene = *buildScene(json, camera);
        
        return new Frame(width, height, raysPerPixel, filePath, backgroundColor, threadCount);
    }

    static FrameServer* buildFrameServer(Camera& camera, Scene& scene, const JSON& json = parseJSON("Settings.json")) {
        int width, height;
        short raysPerPixel, port = 1990;
        Vector backgroundColor = BLACK;
        
        Types* temp = json.get("frame_width");
        if (temp != nullptr) 
            width = (int) temp -> dblValue;

        temp = json.get("frame_height");
        if (temp != nullptr) 
            height = (int) temp -> dblValue;

        temp = json.get("rays_per_pixel");
        if (temp != nullptr) 
            raysPerPixel = (int) temp -> dblValue;

        temp = json.get("server_port");
        if (temp != nullptr) 
            port = (short) temp -> dblValue;

        temp = json.get("background_color");
        if (temp != nullptr)
            backgroundColor = Color(*(temp -> vecValue));

        camera = *buildCamera(json);
        scene = *buildScene(json, camera);
        return new FrameServer(width, height, raysPerPixel, port, backgroundColor);
    }

};

#endif
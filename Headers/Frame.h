#ifndef FRAME
#define FRAME

// Include -lwsock32 when compiling with G++
#include "Vector.h"
#include "Camera.h"
#include "Scene.h"
#include "Render.h"

#include <iostream>
#include <fstream>
#include <WinSock2.h>
#include <string>
#include <filesystem>
#include <thread>
#include <future>
#include <atomic>
#include <mutex>

namespace fs = std::filesystem;
const short colorMult = 255;

class Frame {
    private:
    int width, height, raysPerPixel, threadCount = 10;
    double aliasingMult, aspectRatio;
    Vector backgroundColor;
    FILE* file;

    public:
    Frame(const int& x, const int& y, const int& aliasing, const std::string& fileName, const Vector& bgColor = BLACK, int threadCount = 0) {
        this -> width = x;
        this -> height = y;
        this -> raysPerPixel = aliasing;
        this -> aliasingMult = 1.0 / (double) raysPerPixel;
        this -> threadCount = threadCount;
        this -> backgroundColor = bgColor;
        aspectRatio = (double) width / height;
        file = fopen(fileName.c_str(), "w");

        fprintf(file, "P3\n%d %d\n255\n", width, height);
    }

    int getWidth() {
        return width;
    }  

    int getHeight() {
        return height;
    }

    Vector getBackground() {
        return backgroundColor;
    }

    double getAspectRatio() {
        return aspectRatio;
    }

    int getThreadCount() {
        return threadCount;
    }

    int getAliasing() {
        return raysPerPixel;
    }

    void render(Camera* camera, Scene* scene) {
        const double recWidth = 1.0 / (width - 1), recHeight = 1.0 / (height - 1);
        double offsets[raysPerPixel];

        for (short i = 0; i < raysPerPixel; ++i) {
            offsets[i] = (double) i * aliasingMult;
        }

        if (threadCount > 1) {
            std::vector<std::future<void>> threads;
            threads.reserve(threadCount);

            if (height % threadCount == 0) {
                volatile std::atomic<int> yLimit(0);
                std::vector<std::string> pixelBuffer(threadCount);
                std::mutex threadLock;

                for (int i = 0; i < threadCount; ++i) {
                    threads.emplace_back(std::async([=, &yLimit, &camera, &scene, &offsets, &pixelBuffer, &threadLock]() {
                        std::string res = "";
                        yLimit += height / threadCount; 
                        const int min = yLimit - height / threadCount, max = yLimit;
                        Ray r;
                        Vector sum;
                        
                        for (int y = min; y < max; ++y) {
                            for (int x = 0; x < width; ++x) {
                                sum = Vector(0, 0, 0);

                                for (int a = 0; a < raysPerPixel; ++a) {
                                    r = camera -> generateRay((x + offsets[a]) * recWidth, (y + offsets[a]) * recHeight);
                                    sum += traceRay(r, *scene, backgroundColor).replaceNAN();      
                                }

                                res += (clamp((sum * aliasingMult).sqrt()) * colorMult).str() + "\n";
                            }
                        }

                        pixelBuffer[min * threadCount / height] = res;
                    }));
                }

                for (short currentThread = 0; currentThread < threads.size(); ++currentThread) threads[currentThread].wait(); 
                for (int i = 0; i < threadCount; ++i) {
                    fprintf(file, pixelBuffer[i].c_str());
                    fflush(file);
                }
            } else {
                volatile std::atomic<int> count(0);
                int max = height * width;
                std::vector<Vector> pixelBuffer(max);
                double recYFrame = 1.0 / height;

                for (int i = 0; i < threadCount; ++i) {
                    threads.emplace_back(std::async([=, &count, &camera, &scene, &offsets, &pixelBuffer]() {
                        std::string res = "";
                        int index = 0;
                        Ray r;
                        
                        while (index < max) {
                            index = count++;

                            for (int a = 0; a < raysPerPixel; ++a) {
                                r = camera -> generateRay((index % width + offsets[a]) * recWidth, (index * recYFrame + offsets[a]) * recHeight);
                                pixelBuffer[index] += traceRay(r, *scene, backgroundColor).replaceNAN();      
                            }
                        }
                    })); 
                }

                for (int i = 0; i < threadCount; ++i) threads[i].wait();
                for (int i = 0; i < max; ++i) {
                    write(pixelBuffer[i]);
                }
            }
        } else {
            Vector sum;
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    sum = Vector(0, 0, 0);
                    for (int a = 0; a < raysPerPixel; ++a) {
                        sum += traceRay(camera -> generateRay((x + offsets[a]) * recWidth, (y + offsets[a]) * recHeight), 
                                        *scene, backgroundColor).replaceNAN();
                    }

                    write(sum);
                }
            }
        }

        close();
    }

    void write(const Vector& color) {
        Vector adjustedColor = clamp((color * aliasingMult).sqrt()) * colorMult;

        fprintf(file, "%lf %lf %lf\n", adjustedColor[0], 
                                       adjustedColor[1], 
                                       adjustedColor[2]);
        fflush(file);
    }

    void open(const std::string& fileName) {
        if (file != nullptr)
            fclose(file);
        
        file = fopen(fileName.c_str(), "w");
    }

    void close() {
        _pclose(file);
    }

};

class FrameServer {
    private:
    const char* headers = "HTTP/1.1 200 OK\r\nAccess-Control-Allow-Origin: *\r\nContent-Type: text/plain\r\nConnection: Closed\r\n\r\n";
    int width, height, raysPerPixel;
    double aliasingMult, aspectRatio;
    Vector backgroundColor;

    WSADATA startWSA() {
        WSADATA w;
        int startupResult = WSAStartup(MAKEWORD(2, 2), &w);
        if (startupResult != 0) 
            return startWSA();
        
        return w;
    }

    void bindSock(SOCKET &sock, sockaddr_in &server) {
        int res = bind(sock, (struct sockaddr *) &server, sizeof(server));
        if (res == SOCKET_ERROR) 
            bindSock(sock, server);
    }

    SOCKET createSock(int type) {
        SOCKET sock = socket(type, SOCK_STREAM, IPPROTO_TCP);
        if (sock == INVALID_SOCKET) 
            return createSock(type);
        
        return sock;
    }

    static void shellCommand(std::string cmd, std::string redirect) {
        STARTUPINFOA si = {0};
        PROCESS_INFORMATION pi = {0};
        char* fullCMD;

        if (CreateProcessA(NULL, ("powershell " + cmd + " > " + redirect).data(), NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, NULL, &si, &pi)) {
            WaitForSingleObject(pi.hProcess, INFINITE);
            CloseHandle(pi.hProcess);
            CloseHandle(pi.hThread);
        } else {
            std::cout << GetLastError() << " Shell Command\n";
        }
    }

    static void readCommand(std::string fileRedirect, std::string &result) {
        std::ifstream input(fileRedirect);
        result = "";
        
        while (input.good()) {
            int curr = input.get();
            if (curr >= 32 && curr <= 126) 
                result += (char) curr;
        }
    }

    protected:
    WSADATA wsa;
    SOCKET main, connection;
    struct sockaddr_in info;

    static std::string getIPv4() {
        std::string file;
        shellCommand("ipconfig", "ipaddr.txt");
        readCommand("ipaddr.txt", file);
        std::remove("ipaddr.txt");

        int start = file.find("IPv4");
        while (start < file.size() && file[start - 2] != ':') {
            ++start;
        }

        std::string ip = "";
        do {
            ip += file[start++];
        } while (file[start] == '.' || (file[start] >= '0' && file[start] <= '9'));
        
        std::cout << ip << "\n";
        return ip;
    }

    public:
    FrameServer(const int& width, const int& height, const short& raysPerPixel, const short &port, const Vector& bgColor = BLACK) {
        this -> width = width;
        this -> height = height;
        this -> raysPerPixel = raysPerPixel;
        this -> backgroundColor = bgColor;
        aliasingMult = 1 / raysPerPixel;
        aspectRatio = (double) width / height;

        wsa = startWSA();
        info.sin_family = AF_INET;
        info.sin_addr.S_un.S_addr = inet_addr(getIPv4().c_str());
        info.sin_port = htons(port);

        main = createSock(info.sin_family);
        bindSock(main, info);
    }

    void render(Camera* camera, Scene* scene) {
        const double recWidth = 1.0 / (width - 1), recHeight = 1.0 / (height - 1);
        char buffer[1024];
        listen(main, 1);
        
        while (1) {
            std::cout << "Connection Waiting!\n";
            connection = accept(main, NULL, NULL);
            std::cout << "Connection Accepted!\n";
            recv(connection, buffer, sizeof(buffer), 0);

            std::string message = headers + std::to_string(width) + " " + std::to_string(height) + " " + std::to_string(raysPerPixel);    
            send(connection, message.c_str(), message.length(), 0);
            closesocket(connection);
            Vector color;

            for (short aliasing = 1; aliasing <= raysPerPixel; ++aliasing) {
                message = headers;

                for (int y = 0; y < height; ++y) {
                    for (int x = 0; x < width; ++x) {
                        color = traceRay(camera -> generateRay((x + randDbl()) * recWidth, (y + randDbl()) * recHeight), *scene, backgroundColor);
                        message += color.replaceNAN().str() + "\n";
                    }
                }
                
                connection = accept(main, NULL, NULL);
                recv(connection, buffer, sizeof(buffer), 0);
                send(connection, message.c_str(), message.length(), 0);
                closesocket(connection);
                if ((aliasing / raysPerPixel * 10) % 1 == 0) std::cout << aliasing << " Done!\n";
            }

            return;
        }
    }

    int getWidth() {
        return width;
    }  

    int getHeight() {
        return height;
    }

    Vector getBackground() {
        return backgroundColor;
    }

    double getAspectRatio() {
        return aspectRatio;
    }

    int getAliasing() {
        return raysPerPixel;
    }
    
};

#endif
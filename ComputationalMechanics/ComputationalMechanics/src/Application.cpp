#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <GLFW/glfw3native.h>
#define WINVER 0x0500
#include <windows.h>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include "Renderer.h"

#include "Module.h"

#include "Problem1.h"
#include "Problem2.h"

module::Module* currentModule = nullptr;

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    ImGuiIO& io = ImGui::GetIO();
    if (!io.WantCaptureMouse)
    {
        currentModule->KeyCallback(window, key, scancode, action, mods);
    }
}

void CursorPosCallback(GLFWwindow* window, double xpos, double ypos)
{
    ImGuiIO& io = ImGui::GetIO();
    if (!io.WantCaptureMouse)
    {
        currentModule->CursorPosCallback(window, xpos, ypos);
    }
}

void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) 
{
    ImGuiIO& io = ImGui::GetIO();
    if (!io.WantCaptureMouse)
    {
        currentModule->MouseButtonCallback(window, button, action, mods);
    }
}

void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset) 
{
    
    ImGuiIO& io = ImGui::GetIO();
    if (!io.WantCaptureMouse)
    {
        currentModule->ScrollCallback(window, xoffset, yoffset);
    }
}

int main(void)
{
    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(800, 600, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, KeyCallback);
    glfwSetCursorPosCallback(window, CursorPosCallback);
    glfwSetMouseButtonCallback(window, MouseButtonCallback);
    glfwSetScrollCallback(window, ScrollCallback);

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    glfwSwapInterval(1);//vertical synchronization

    if(glewInit() != GLEW_OK)
        std::cout<<"glewinit error"<<std::endl;
    std::cout << glfwGetVersionString() << std::endl;

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    {//make a scope to clean the vartables in stack,otherwise somethings go wrong after call glfwTerminate();

        module::ModuleMenu* moduleMenu = new module::ModuleMenu(currentModule);
        currentModule = moduleMenu;

        double deltaTime = 0.0f;
        double lastFrame = 0.0f;

        moduleMenu->RegisterModule<module::Problem1>("Problem1");
        moduleMenu->RegisterModule<module::Problem2>("Problem2");

        /* Loop until the user closes the window */
        while (!glfwWindowShouldClose(window))
        {
            double currentFrame = glfwGetTime();
            deltaTime = currentFrame - lastFrame;
            lastFrame = currentFrame;

            currentModule->OnUpdate(deltaTime);

            currentModule->OnRender();

            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            {
                ImGui::Begin("test");

                if (currentModule != moduleMenu && ImGui::Button("<-"))
                {
                    delete currentModule;
                    currentModule = moduleMenu;
                }
                currentModule->OnImguiRender();

                ImGui::End();
            }
            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

            /* Swap front and back buffers */
            glfwSwapBuffers(window);

            /* Poll for and process events */
            glfwPollEvents();
        }
        if (currentModule != moduleMenu)
            delete currentModule;
        delete moduleMenu;
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwTerminate();
    return 0;
}
#include "Problem2.h"

#include "Renderer.h"
#include "VertexBuffer.h"
#include "IndexBuffer.h"
#include "VertexArray.h"

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <fstream>

static const double pi = 3.141592653589793;

namespace module {

    Problem2::Problem2()
        :
        m_Proj(glm::mat4(1.0f)),
        m_View(glm::mat4(1.0f)),
        m_Model(glm::mat4(1.0f)),
        m_MVP(glm::mat4(1.0f)),
        m_IO(ImGui::GetIO()),
        m_scale(0.0f),
        n1(16),
        n2(16),
        insideBoundary(BoundaryType::FixedTemperature),
        insideBoundaryValue(400),
        outsideBoundary(BoundaryType::FixedTemperature),
        outsideBoundaryValue(-100)
    {
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);

        m_VAO = std::make_unique<VertexArray>();
        m_VBO = std::make_unique<VertexBuffer>(m_Vertices.data(), sizeof(Vertex) * m_Vertices.size());
        VertexBufferLayout layout;
        layout.Push<float>(3);//position
        layout.Push<float>(3);//color
        layout.Vacate(sizeof(double));//temperature
        m_VAO->AddBuffer(*m_VBO, layout);
        m_IBO = std::make_unique<IndexBuffer>(m_Indices.data(), m_Indices.size());
        m_Shader = std::make_unique<Shader>("res/shaders/Problem2.shader");
        m_Shader->Bind();
        m_Camera = std::make_unique<Camera>(m_View);

        glEnable(GL_POLYGON_OFFSET_FILL); // 启用填充面的多边形偏移
        glPolygonOffset(0.5f, 0.1f); // 设置多边形偏移因子和单位

        m_Camera->m_cameraPos = { -0.58f, 1.42f, 0.41f };
        m_Camera->yaw = 250.0f;
        m_Camera->pitch = -48.0f;
        GenerateMesh();
    }

    Problem2::~Problem2()
    {
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_POLYGON_OFFSET_FILL);
    }

    void Problem2::OnUpdate(double deltaTime)
    {
        m_Camera->CameraUpdate(deltaTime);
        m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(pow(10.0f, m_scale)));//调整模型大小
    }

    void Problem2::OnRender()
    {
        int width, height;
        GLFWwindow* window = glfwGetCurrentContext();
        if (!glfwGetWindowAttrib(window, GLFW_ICONIFIED))
        {
            //窗口没有最小化
            glfwGetWindowSize(window, &width, &height);
            m_Proj = glm::perspective(glm::radians(m_Camera->fov), (float)width / (float)height, 0.001f, 1000.0f);

        }
        m_MVP = m_Proj * m_View * m_Model;
        m_Shader->SetUniformMat4f("u_MVP", m_MVP);

        /* Render here */
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        Renderer renderer;
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        m_VBO->ReData(m_Vertices.data(), sizeof(Vertex) * m_Vertices.size());
        m_IBO->ReData(m_Indices.data(), m_Indices.size());

        glLineWidth(1.0f);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        m_Shader->SetUniform1i("u_Mode", 1);
        renderer.DrawTriangle(*m_VAO, *m_IBO, *m_Shader, m_Indices.size());
        m_Shader->SetUniform1i("u_Mode", 0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        renderer.DrawTriangle(*m_VAO, *m_IBO, *m_Shader, m_Indices.size());
    }

    void Problem2::OnImguiRender()
    {
        ImGui::SliderInt("mesh size", &n1, 8, 96);

        if (ImGui::Button("Resize Mesh"))
        {
            GenerateMesh();
        }
        if (ImGui::Button("Calculate"))
        {
            TemperatureUpdate();
            ColorUpdate();
        }
        if (ImGui::Button("PrintTemperature"))
        {
            PrintTemperature();
        }
        if (ImGui::Button("OutputData"))
        {
            OutputData("Temperature.data");
        }

        ImGui::Text("insideBoundary:");
        ImGui::PushItemWidth(50);  // 将输入框的宽度设置为 100 像素
        ImGui::InputText("value1", m_input1, IM_ARRAYSIZE(m_input1));
        ImGui::PopItemWidth();  // 恢复默认的宽度
        ImGui::SameLine();
        ImGui::PushID("insideBoundary");//防止ID冲突
        ImGui::RadioButton("FixedTemperature",  (int*)&insideBoundary, 0);
        ImGui::SameLine();
        ImGui::RadioButton("FixedHeatFlux", (int*)&insideBoundary, 1);
        ImGui::PopID();
        

        ImGui::Text("outsideBoundary:");
        ImGui::PushItemWidth(50);
        ImGui::InputText("value2", m_input2, IM_ARRAYSIZE(m_input2));
        ImGui::PopItemWidth();
        ImGui::SameLine();
        ImGui::PushID("outsideBoundary");
        ImGui::RadioButton("FixedTemperature", (int*)&outsideBoundary, 0);
        ImGui::SameLine();
        ImGui::RadioButton("FixedHeatFlux", (int*)&outsideBoundary, 1);
        ImGui::PopID();

        insideBoundaryValue = atof(m_input1);
        outsideBoundaryValue = atof(m_input2);
        n2 = n1;

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / m_IO.Framerate, m_IO.Framerate);

    }

    void Problem2::CursorPosCallback(GLFWwindow* window, double xpos, double ypos)
    {
        Camera::CursorPosCallback(window, xpos, ypos);
    }

    void Problem2::KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
    {
        Camera::KeyCallback(window, key, scancode, action, mods);
    }

    void Problem2::MouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
    {
        Camera::MouseButtonCallback(window, button, action, mods);
    }

    void Problem2::ScrollCallback(GLFWwindow* window, double xoffset, double yoffset)
    {
        Camera::ScrollCallback(window, xoffset, yoffset);
    }

    void Problem2::GenerateMesh()
    {
        //生成n1*n2*2个三角网格
        m_Vertices.clear();
        m_Faces.clear();
        m_Indices.clear();
        double ra = 0.5;
        double rb = 1.5;
        //创建顶点
        for (int i = 0;i < n1 + 1;i++)
        {
            for (int j = 0;j < n2 + 1;j++)
            {
                double rho = ra + i * (rb - ra) / n1;
                double theta = j * pi / 2 / n2;
                double x = rho * std::cos(theta);
                double y = rho * std::sin(theta);
                m_Vertices.push_back({ { x,y,0 } , {1.0f, 1.0f, 1.0f}, 0.0 });
            }
        }
        //创建三角形
        for (int i = 0;i < n1;i++)
        {
            for (int j = 0;j < n2;j++)
            {
                Face f1;
                Face f2;
                f1.v1 = (n2 + 1) * i + j;
                f1.v2 = (n2 + 1) * (i + 1) + j;
                f1.v3 = (n2 + 1) * (i + 1) + j + 1;
                f2.v1 = (n2 + 1) * i + j;
                f2.v2 = (n2 + 1) * (i + 1) + j + 1;
                f2.v3 = (n2 + 1) * i + j + 1;
                m_Indices.push_back(f1.v1);
                m_Indices.push_back(f1.v2);
                m_Indices.push_back(f1.v3);
                m_Indices.push_back(f2.v1);
                m_Indices.push_back(f2.v2);
                m_Indices.push_back(f2.v3);
                m_Faces.push_back(f1);
                m_Faces.push_back(f2);
            }
        }
        //计算每个三角形的矩阵K
        for (auto& face : m_Faces)
        {
            glm::vec3 v1 = m_Vertices[face.v1].position;
            glm::vec3 v2 = m_Vertices[face.v2].position;
            glm::vec3 v3 = m_Vertices[face.v3].position;
            double delta = (-(v1.x - v2.x) * (v3.y - v2.y) + (v1.y - v2.y) * (v3.x - v2.x)) / 2;
            double a1 = v2.y - v3.y;
            double b1 = v3.x - v2.x;
            double a2 = v3.y - v1.y;
            double b2 = v1.x - v3.x;
            double a3 = v1.y - v2.y;
            double b3 = v2.x - v1.x;
            double k = 1;//热传导系数
            Eigen::MatrixXd matrixB(2, 3);
            matrixB << a1, a2, a3,
                b1, b2, b3;
            face.matrixK = k / 4 / delta * matrixB.transpose() * matrixB;
        }
    }
    void Problem2::ClearTemperature(double initialTemperature)
    {
        for (auto& vertex : m_Vertices)
        {
            vertex.temperature = initialTemperature;
        }
    }

    void Problem2::TemperatureUpdate()
    {
        Eigen::SparseMatrix<double> matrixA((n1 + 1) * (n2 + 1), (n1 + 1) * (n2 + 1));
        Eigen::MatrixXd matrixb((n1 + 1) * (n2 + 1), 1);
        matrixb.setZero();
        for (auto& face : m_Faces)
        {
            matrixA.coeffRef(face.v1, face.v1) += face.matrixK.coeff(0, 0);
            matrixA.coeffRef(face.v1, face.v2) += face.matrixK.coeff(0, 1);
            matrixA.coeffRef(face.v1, face.v3) += face.matrixK.coeff(0, 2);
            matrixA.coeffRef(face.v2, face.v1) += face.matrixK.coeff(1, 0);
            matrixA.coeffRef(face.v2, face.v2) += face.matrixK.coeff(1, 1);
            matrixA.coeffRef(face.v2, face.v3) += face.matrixK.coeff(1, 2);
            matrixA.coeffRef(face.v3, face.v1) += face.matrixK.coeff(2, 0);
            matrixA.coeffRef(face.v3, face.v2) += face.matrixK.coeff(2, 1);
            matrixA.coeffRef(face.v3, face.v3) += face.matrixK.coeff(2, 2);
        }
        //处理边界条件，如果是热流边界条件，则修改matrixb，相当于载荷，如果是温度边界条件，则矩阵降阶（或修改矩阵A和b），相当于固定顶点
        //内璧边界条件
        if (insideBoundary == BoundaryType::FixedTemperature)
        {
            for (int i = 0;i < n2 + 1;i++)
            {
                for (int j = 0;j < (n1 + 1) * (n2 + 1);j++)
                {
                    matrixb.coeffRef(j, 0) -= matrixA.coeff(j, i) * insideBoundaryValue;
                    matrixA.coeffRef(i, j) = 0;
                    matrixA.coeffRef(j, i) = 0;
                }
                matrixA.coeffRef(i, i) = 1;
                matrixb.coeffRef(i, 0) = insideBoundaryValue;
            }
        }
        else
        {
            for (int i = 1;i < n2;i++)
            {
                matrixb.coeffRef(i, 0) += insideBoundaryValue * (0.5 * pi * 0.5) / n2;
            }
            matrixb.coeffRef(0, 0) += insideBoundaryValue * (0.5 * pi * 0.5) / n2 / 2;
            matrixb.coeffRef(n2, 0) += insideBoundaryValue * (0.5 * pi * 0.5) / n2 / 2;
        }
        //外璧边界条件
        if (outsideBoundary == BoundaryType::FixedTemperature)
        {
            for (int i = n1 * (n2 + 1);i < (n1 + 1) * (n2 + 1);i++)
            {
                for (int j = 0;j < (n1 + 1) * (n2 + 1);j++)
                {
                    matrixb.coeffRef(j, 0) -= matrixA.coeff(j, i) * outsideBoundaryValue;
                    matrixA.coeffRef(i, j) = 0;
                    matrixA.coeffRef(j, i) = 0;
                }
                matrixA.coeffRef(i, i) = 1;
                matrixb.coeffRef(i, 0) = outsideBoundaryValue;
            }
        }
        else
        {
            for (int i = n1 * (n2 + 1) + 1;i < (n1 + 1) * (n2 + 1) - 1;i++)
            {
                matrixb.coeffRef(i, 0) += outsideBoundaryValue * (0.5 * pi * 1.5) / n2;
            }
            matrixb.coeffRef(n1 * (n2 + 1), 0) += outsideBoundaryValue * (0.5 * pi * 1.5) / n2 / 2;
            matrixb.coeffRef((n1 + 1) * (n2 + 1) - 1, 0) += outsideBoundaryValue * (0.5 * pi * 1.5) / n2 / 2;
        }
        //侧面边界条件
        //相当于热流边界条件，但是值为0，无需操作

        //解方程Ax=b
        //是方阵，能用LU分解
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(matrixA);
        solver.factorize(matrixA);
        Eigen::MatrixXd matrixx((n1 + 1) * (n2 + 1), 1);
        matrixx = solver.solve(matrixb);

        for (int i = 0;i < (n1 + 1) * (n2 + 1);i++)
        {
            m_Vertices[i].temperature = matrixx.coeff(i, 0);
        }
    }

    void Problem2::ColorUpdate()
    {
        for (auto& vertex : m_Vertices)
        {
            glm::vec3 colorMin = glm::vec3(1.0f, 1.0f, 0.0f);//黄色
            glm::vec3 colorMax = glm::vec3(1.0f, 0.0f, 0.0f);//红色
            double valueMin = 300.0;
            double valueMax = 400.0;
            float ratio;
            if (vertex.temperature <= valueMin)
            {
                ratio = 0;
            }
            else if (vertex.temperature >= valueMax)
            {
                ratio = 1;
            }
            else
            {
                ratio = (vertex.temperature - valueMin) / (valueMax - valueMin);
            }
            vertex.color = ratio * colorMax + (1 - ratio) * colorMin;
        }
    }

    void Problem2::PrintTemperature()
    {
        for (int i = 0;i < n2 + 1;i++)
        {
            for (int j = 0;j < n1 + 1;j++)
            {
                printf("v%d(%d,%d):\t%.3f\n", (n1 + 1) * i + j, i, j, m_Vertices[(n1 + 1) * i + j].temperature);
            }
        }
    }

    void Problem2::PrintMatrixK()
    {
        for (int i = 0;i < m_Faces.size();i++)
        {
            std::cout << "face:" << i << std::endl << m_Faces[i].matrixK << std::endl;
        }
    }

    void Problem2::OutputData(std::string filepath)
    {
        //产生用于绘图的数据，输出到外部文件里
        std::ofstream fout;
        fout.open(filepath, std::ios::out);
        float x;
        float T;
        for (int i = 0;i < n1 + 1;i++)
        {
            x = 0.5 * (float)i / n1;
            T = m_Vertices[(n2 + 1) * i].temperature;
            fout << x << " " << T << std::endl;
        }
        fout.close();

        std::cout << "output data" << std::endl;
    }
}
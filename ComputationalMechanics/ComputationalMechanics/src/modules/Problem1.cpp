#include "Problem1.h"

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


namespace module {

    Problem1::Problem1()
        :
        m_Proj(glm::mat4(1.0f)),
        m_View(glm::mat4(1.0f)),
        m_Model(glm::mat4(1.0f)),
        m_MVP(glm::mat4(1.0f)),
        m_IO(ImGui::GetIO()),
        m_WireframeMode(false),
        m_scale(-2.0f),
        /*E(100 * 4.45 / pow(0.0254, 2)),
		A(0.2 * pow(0.0254, 2)),
		G(38.5 * 4.45 / pow(0.0254, 2)),
		J(0.006 * pow(0.0254, 4)),
		Iy(0.003 * pow(0.0254, 4)),
		Iz(0.003 * pow(0.0254, 4)),*/
        E(100),
        A(0.2),
        G(38.5),
        J(0.006),
        Iy(0.003),
        Iz(0.003)
    {
        
        std::vector<glm::vec3> Vertices =
        {
            { 58.0f, 38.0f, 0.0f  },//0  1
            { 48.0f, 38.0f, 0.0f  },//1  2
            { 31.0f, 38.0f, 0.0f  },//2  3
            { 17.0f, 38.0f, 22.0f },//3  4
            { 0.0f,  38.0f, 24.0f },//4  5
            { 58.0f, 38.0f, 42.0f },//5  6
            { 48.0f, 38.0f, 42.0f },//6  7
            { 36.0f, 38.0f, 70.0f },//7  8
            { 0.0f,  38.0f, 75.0f },//8  9

            { 58.0f, -38.0f, 0.0f  },//9   1'
            { 48.0f, -38.0f, 0.0f  },//10  2'
            { 31.0f, -38.0f, 0.0f  },//11  3'
            { 17.0f, -38.0f, 22.0f },//12  4'
            { 0.0f,  -38.0f, 24.0f },//13  5'
            { 58.0f, -38.0f, 42.0f },//14  6'
            { 48.0f, -38.0f, 42.0f },//15  7'
            { 36.0f, -38.0f, 70.0f },//16  8'
            { 0.0f,  -38.0f, 75.0f },//17  9'

            { 58.0f, 17.0f, 42.0f },//18  10
            { 58.0f, 17.0f, 0.0f  },//19  11
            { 0.0f,  17.0f, 0.0f  },//20  12
            { 0.0f,  17.0f, 24.0f },//21  13

            { 58.0f, -17.0f, 42.0f },//22  10'
            { 58.0f, -17.0f, 0.0f  }, //23  11'
            { 0.0f,  -17.0f, 0.0f  }, //24  12'
            { 0.0f,  -17.0f, 24.0f },//25  13'

            { 18.0f, 0.0f,  72.0f },//26  14
            { 0.0f,  0.0f,  37.5f },//27  15
        };
        
        /*std::vector<glm::vec3> Vertices =
        {
            {10.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f},
        };*/
        m_Vertices.insert(m_Vertices.end(), Vertices.begin(), Vertices.end());
        
        std::vector<unsigned int> Indices =
        {
            0, 1,
            1, 2,
            2, 3,
            3, 4,
            4, 8,
            5, 6,
            6, 7,
            7, 8,
            1, 6,

            9, 10,
            10, 11,
            11, 12,
            12, 13,
            13, 17,
            14, 15,
            15, 16,
            16, 17,
            10, 15,

            0, 19,
            5, 18,
            9, 23,
            14, 22,
            19, 23,
            18, 22,
            18, 19,
            19, 20,
            20, 21,
            22, 23,
            23, 24,
            24, 25,


            26, 7,
            26, 8,
            26, 16,
            26, 17,
            27, 8,
            27, 4,
            27, 17,
            27, 13,
            7, 16,
            8, 17,
            4, 21,
            21, 25,
            25, 13,
        };
        
        /*std::vector<unsigned int> Indices =
        {
            0,1,
        };*/
        m_Indices.insert(m_Indices.end(), Indices.begin(), Indices.end());

        //���ĸ����ǹ̶��ģ�����ά����������λ��
        m_FixedVertices.push_back(19);
        m_FixedVertices.push_back(20);
        m_FixedVertices.push_back(23);
        m_FixedVertices.push_back(24);
        //m_FixedVertices.push_back(1);

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);

        m_VAO = std::make_unique<VertexArray>();
        m_VBO = std::make_unique<VertexBuffer>(m_Vertices.data(), sizeof(glm::vec3) * m_Vertices.size());
        VertexBufferLayout layout;
        layout.Push<float>(3);//position
        m_VAO->AddBuffer(*m_VBO, layout);
        m_IBO = std::make_unique<IndexBuffer>(m_Indices.data(), m_Indices.size());
        m_Shader = std::make_unique<Shader>("res/shaders/Lines.shader");
        m_Shader->Bind();
        m_Camera = std::make_unique<Camera>(m_View);

        //glEnable(GL_POLYGON_OFFSET_FILL); // ���������Ķ����ƫ��
        //glPolygonOffset(0.5f, 0.1f); // ���ö����ƫ�����Ӻ͵�λ
        printf("num:%d\n", m_Indices.size());
    }

    Problem1::~Problem1()
    {
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_POLYGON_OFFSET_FILL);
    }

    void Problem1::OnUpdate(double deltaTime)
    {
        m_Camera->CameraUpdate(deltaTime);
        m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(pow(10.0f, m_scale)));//����ģ�ʹ�С
    }

    void Problem1::OnRender()
    {
        if (m_WireframeMode)
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        else
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        int width, height;
        GLFWwindow* window = glfwGetCurrentContext();
        if (!glfwGetWindowAttrib(window, GLFW_ICONIFIED))
        {
            //����û����С��
            glfwGetWindowSize(window, &width, &height);
            m_Proj = glm::perspective(glm::radians(m_Camera->fov), (float)width / (float)height, 0.1f, 1000.0f);

        }
        m_MVP = m_Proj * m_View * m_Model;
        m_Shader->SetUniformMat4f("u_MVP", m_MVP);

        /* Render here */
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        Renderer renderer;
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        m_VBO->ReData(m_Vertices.data(), sizeof(glm::vec3) * m_Vertices.size());
        m_IBO->ReData(m_Indices.data(), m_Indices.size());

        glLineWidth(5.0f);
        m_Shader->Bind();
        m_VAO->Bind();
        m_VBO->Bind();
        glDrawElements(GL_LINES, m_Indices.size(), GL_UNSIGNED_INT, nullptr);
    }

    void Problem1::OnImguiRender()
    {
        if (ImGui::Button("Calculate"))
        {
            Calculate();
        }

        ImGui::SliderFloat("model scale", &m_scale, -2.0f, 2.0f);

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / m_IO.Framerate, m_IO.Framerate);

    }

    void Problem1::CursorPosCallback(GLFWwindow* window, double xpos, double ypos)
    {
        Camera::CursorPosCallback(window, xpos, ypos);
    }

    void Problem1::KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
    {
        Camera::KeyCallback(window, key, scancode, action, mods);
    }

    void Problem1::MouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
    {
        Camera::MouseButtonCallback(window, button, action, mods);
    }

    void Problem1::ScrollCallback(GLFWwindow* window, double xoffset, double yoffset)
    {
        Camera::ScrollCallback(window, xoffset, yoffset);
    }

    void Problem1::Calculate()
    {
        // Ax = b
        Eigen::SparseMatrix<double> matrixA(m_Vertices.size() * 6, m_Vertices.size() * 6);
        Eigen::MatrixXd matrixb(m_Vertices.size() * 6, 1);
        for (int i = 0;i < m_Indices.size() / 2;i++)
        {
            //��ÿ����ѭ��
            int v1Index = m_Indices[2 * i];
            int v2Index = m_Indices[2 * i + 1];
            glm::vec3 v1 = m_Vertices[v1Index];
            glm::vec3 v2 = m_Vertices[v2Index];
            glm::vec3 diretion = v2 - v1;
            glm::vec3 axis_x;//�ֲ�����ϵ����������������ϵ�е�����
            glm::vec3 axis_y;
            glm::vec3 axis_z;
            double l = glm::length(diretion);
            if (diretion.x == 0 && diretion.y == 0)
            {
                //����ֱ�������⴦��
                axis_x = glm::vec3(0.0f, 0.0f, 1.0f);
                axis_y = glm::vec3(1.0f, 0.0f, 0.0f);
                axis_z = glm::vec3(0.0f, 1.0f, 0.0f);
            }
            else if (diretion.z == 0)
            {
                //��ˮƽ�������⴦��
                axis_x = glm::normalize(diretion);
                axis_z = glm::vec3(0.0f, 0.0f, 1.0f);
                axis_y = glm::normalize(glm::cross(axis_z, axis_x));
            }
            else
            {
                axis_x = glm::normalize(diretion);
                glm::vec3 projection = glm::vec3(axis_x.x, axis_x.y, 0.0f);//diretion��ˮƽ��ͶӰ
                axis_y = glm::normalize(glm::cross(axis_x, projection));
                axis_z = glm::normalize(glm::cross(axis_x, axis_y));
            }
            Eigen::MatrixXd matrix_lambda(3, 3);//����ƴ������任����
            {
                glm::vec3 x = glm::vec3(1.0f, 0.0f, 0.0f);
                glm::vec3 y = glm::vec3(0.0f, 1.0f, 0.0f);
                glm::vec3 z = glm::vec3(0.0f, 0.0f, 1.0f);
                matrix_lambda <<
                    (double)glm::dot(x, axis_x), (double)glm::dot(x, axis_y), (double)glm::dot(x, axis_z),
                    (double)glm::dot(y, axis_x), (double)glm::dot(y, axis_y), (double)glm::dot(y, axis_z),
                    (double)glm::dot(z, axis_x), (double)glm::dot(z, axis_y), (double)glm::dot(z, axis_z);
            }
            Eigen::MatrixXd matrixT(12, 12);//����任����
            matrixT.setZero();
            matrixT.block(0, 0, matrix_lambda.rows(), matrix_lambda.cols()) = matrix_lambda;
            matrixT.block(3, 3, matrix_lambda.rows(), matrix_lambda.cols()) = matrix_lambda;
            matrixT.block(6, 6, matrix_lambda.rows(), matrix_lambda.cols()) = matrix_lambda;
            matrixT.block(9, 9, matrix_lambda.rows(), matrix_lambda.cols()) = matrix_lambda;
            Eigen::MatrixXd matrixK_(12, 12);//�ֲ�����ϵ�µĸնȾ���
            Eigen::MatrixXd matrixK(12, 12);//��������ϵ�µĸնȾ���
            //��matrixK_����Ԫ��
            {
                matrixK_.setZero();

                matrixK_.coeffRef(0, 0) = (double)E * A / l;
                matrixK_.coeffRef(0, 6) = (double)-E * A / l;
                matrixK_.coeffRef(6, 0) = (double)-E * A / l;
                matrixK_.coeffRef(6, 6) = (double)E * A / l;

                matrixK_.coeffRef(3, 3) = (double)G * J / l;
                matrixK_.coeffRef(3, 9) = (double)-G * J / l;
                matrixK_.coeffRef(9, 3) = (double)-G * J / l;
                matrixK_.coeffRef(9, 9) = (double)G * J / l;

                matrixK_.coeffRef(1, 1)  = (double)12 * E * Iz / pow(l, 3);
                matrixK_.coeffRef(1, 5)  = (double)6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(1, 7)  = (double)-12 * E * Iz / pow(l, 3);
                matrixK_.coeffRef(1, 11) = (double)6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(5, 1)  = (double)6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(5, 5)  = (double)4 * E * Iz / pow(l, 1);
                matrixK_.coeffRef(5, 7)  = (double)-6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(5, 11) = (double)2 * E * Iz / pow(l, 1);
                matrixK_.coeffRef(7, 1)  = (double)-12 * E * Iz / pow(l, 3);
                matrixK_.coeffRef(7, 5)  = (double)-6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(7, 7)  = (double)12 * E * Iz / pow(l, 3);
                matrixK_.coeffRef(7, 11) = (double)-6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(11, 1) = (double)6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(11, 5) = (double)2 * E * Iz / pow(l, 1);
                matrixK_.coeffRef(11, 7) = (double)-6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(11, 11)= (double)4 * E * Iz / pow(l, 1);

                matrixK_.coeffRef(2, 2)  = (double)-12 * E * Iz / pow(l, 3);
                matrixK_.coeffRef(2, 4)  = (double)6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(2, 8)  = (double)12 * E * Iz / pow(l, 3);
                matrixK_.coeffRef(2, 10) = (double)6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(4, 2)  = (double)6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(4, 4)  = (double)-4 * E * Iz / pow(l, 1);
                matrixK_.coeffRef(4, 8)  = (double)-6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(4, 10) = (double)-2 * E * Iz / pow(l, 1);
                matrixK_.coeffRef(8, 2)  = (double)12 * E * Iz / pow(l, 3);
                matrixK_.coeffRef(8, 4)  = (double)-6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(8, 8)  = (double)-12 * E * Iz / pow(l, 3);
                matrixK_.coeffRef(8, 10) = (double)-6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(10, 2) = (double)6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(10, 4) = (double)-2 * E * Iz / pow(l, 1);
                matrixK_.coeffRef(10, 8) = (double)-6 * E * Iz / pow(l, 2);
                matrixK_.coeffRef(10, 10)= (double)-4 * E * Iz / pow(l, 1);
            }
            //K = T(T) * K_ * T
            matrixK = matrixT.transpose() * matrixK_ * matrixT;
            //��K(12x12)�ֳ�4��(6x6)�ľ����볬��ϡ�����A
            for (int j = 0;j < 6;j++)
            {
                for (int k = 0;k < 6;k++)
                {
                    if (matrixK.coeff(j, k) != 0)
                    {
                        matrixA.coeffRef(v1Index * 6 + j, v1Index * 6 + k) += matrixK.coeff(j, k);
                    }

                }
            }
            for (int j = 0;j < 6;j++)
            {
                for (int k = 0;k < 6;k++)
                {
                    if (matrixK.coeff(j, k + 6) != 0)
                    {
                        matrixA.coeffRef(v1Index * 6 + j, v2Index * 6 + k) += matrixK.coeff(j, k + 6);
                    }

                }
            }
            for (int j = 0;j < 6;j++)
            {
                for (int k = 0;k < 6;k++)
                {
                    if (matrixK.coeff(j + 6, k) != 0)
                    {
                        matrixA.coeffRef(v2Index * 6 + j, v1Index * 6 + k) += matrixK.coeff(j + 6, k);
                    }

                }
            }
            for (int j = 0;j < 6;j++)
            {
                for (int k = 0;k < 6;k++)
                {
                    if (matrixK.coeff(j + 6, k + 6) != 0)
                    {
                        matrixA.coeffRef(v2Index * 6 + j, v2Index * 6 + k) += matrixK.coeff(j + 6, k + 6);
                    }

                }
            }
            if(i == 0)
            std::cout << matrixK << std::endl<<std::endl;
        }//��ÿ������ѭ��

        //����Ŀǰ�غ��������٣���û�иı��������ֱ�������ﹹ�����b����������ѭ�������Ԫ��
        matrixb.setZero();
        matrixb.coeffRef(0, 0) = -3194.0;
        matrixb.coeffRef(1, 0) = -856.0;

        //���ݹ̶���߽������޸ľ���
        for (int j = 0;j < m_FixedVertices.size();j++)
        {
            printf("index = %d\n", m_FixedVertices[j]);
            int vIndex = m_FixedVertices[j];//Ҫ�̶��ĵ������
            for (int k = 0;k < m_Vertices.size() * 6;k++)
            {
                matrixA.coeffRef(vIndex * 6 + 0, k) = 0.0;
                matrixA.coeffRef(vIndex * 6 + 1, k) = 0.0;
                matrixA.coeffRef(vIndex * 6 + 2, k) = 0.0;
                matrixA.coeffRef(vIndex * 6 + 3, k) = 0.0;
                matrixA.coeffRef(vIndex * 6 + 4, k) = 0.0;
                matrixA.coeffRef(vIndex * 6 + 5, k) = 0.0;

                matrixA.coeffRef(k, vIndex * 6 + 0) = 0.0;
                matrixA.coeffRef(k, vIndex * 6 + 1) = 0.0;
                matrixA.coeffRef(k, vIndex * 6 + 2) = 0.0;
                matrixA.coeffRef(k, vIndex * 6 + 3) = 0.0;
                matrixA.coeffRef(k, vIndex * 6 + 4) = 0.0;
                matrixA.coeffRef(k, vIndex * 6 + 5) = 0.0;

            }
            matrixA.coeffRef(vIndex * 6 + 0, vIndex * 6 + 0) = 1.0;
            matrixA.coeffRef(vIndex * 6 + 1, vIndex * 6 + 1) = 1.0;
            matrixA.coeffRef(vIndex * 6 + 2, vIndex * 6 + 2) = 1.0;
            matrixA.coeffRef(vIndex * 6 + 3, vIndex * 6 + 3) = 1.0;
            matrixA.coeffRef(vIndex * 6 + 4, vIndex * 6 + 4) = 1.0;
            matrixA.coeffRef(vIndex * 6 + 5, vIndex * 6 + 5) = 1.0;
            
            matrixb.coeffRef(vIndex * 6 + 0, 0) = 0.0;
            matrixb.coeffRef(vIndex * 6 + 1, 0) = 0.0;
            matrixb.coeffRef(vIndex * 6 + 2, 0) = 0.0;
            matrixb.coeffRef(vIndex * 6 + 3, 0) = 0.0;
            matrixb.coeffRef(vIndex * 6 + 4, 0) = 0.0;
            matrixb.coeffRef(vIndex * 6 + 5, 0) = 0.0;
            //�������ڱ߽�̶�����0������matrixb����������Ҫ��ȥ�̶�ֵ
        }


        //debug
        /*for (int k = 0;k < m_Vertices.size() * 6;k++)
        {
            printf("(120,%d):%.2lf,%.2lf\n", k,matrixA.coeff(120, k), matrixb.coeff(k, 0));
        }*/


        //�ⷽ��Ax=b
        //�Ƿ�������LU�ֽ�
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(matrixA);
        solver.factorize(matrixA);
        matrixx = solver.solve(matrixb);

        /*Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr;
        matrixA.makeCompressed();
        qr.compute(matrixA);
        if (qr.info() != Eigen::Success) 
        {
            std::cout << "Sparse QR decomposition failed." << std::endl;
            return;
        }
        Eigen::MatrixXd matrixx(m_Vertices.size() * 6, 1);
        matrixx = qr.solve(matrixb);*/

        //������
        for (int j = 0;j < m_Vertices.size();j++)
        {
            printf("\nvertex %d:\n",j);
            printf("%.2f,\t\t%.2f,\t\t%.2f\n", matrixx.coeff(j * 6 + 0, 0), matrixx.coeff(j * 6 + 1, 0), matrixx.coeff(j * 6 + 2, 0));
            printf("%.2f,\t\t%.2f,\t\t%.2f\n", matrixx.coeff(j * 6 + 3, 0), matrixx.coeff(j * 6 + 4, 0), matrixx.coeff(j * 6 + 5, 0));
        }
        
    }
}
#pragma once

#include "Module.h"

#include "Shader.h"
#include "VertexBufferLayout.h"
#include "Texture.h"
#include "Camera.h"
#include "HalfEdge2.h"

#include "GLFW/glfw3.h"

#include "imgui/imgui.h"

#include <memory>
#include <unordered_set>

#include "obj-loader/OBJ-Loader.h"

namespace module {

	class Problem1 : public Module
	{
	public:
		Problem1();
		~Problem1();

		void OnUpdate(double deltaTime) override;
		void OnRender() override;
		void OnImguiRender() override;
		void CursorPosCallback(GLFWwindow* window, double xpos, double ypos) override;
		void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) override;
		void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) override;
		void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset) override;

		void CalculateDeflections();
		void CalculateMoments();
	private:

		glm::mat4 m_Proj;
		glm::mat4 m_View;
		glm::mat4 m_Model;
		glm::mat4 m_MVP;

		std::unique_ptr<VertexArray> m_VAO;
		std::unique_ptr<VertexBuffer> m_VBO;
		std::unique_ptr<IndexBuffer> m_IBO;
		std::unique_ptr<Shader> m_Shader;
		std::unique_ptr<Camera> m_Camera;

		std::vector<glm::vec3> m_Vertices;
		std::vector<unsigned int> m_Indices;
		std::vector<int> m_FixedVertices;

		float m_scale;

		double E;
		double A;
		double G;
		double J;
		double Iy;
		double Iz;
		Eigen::MatrixXd matrixx;

		ImGuiIO& m_IO;
	};
}
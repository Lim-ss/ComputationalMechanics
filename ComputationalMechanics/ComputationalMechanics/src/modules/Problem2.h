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

static struct Face
{
	int v1;
	int v2;
	int v3;
	Eigen::Matrix3d matrixK;
};

static struct Vertex
{
	glm::vec3 position;
	glm::vec3 color;
	double temperature;
};

static enum BoundaryType
{
	FixedTemperature = 0,
	FixedHeatFlux = 1,
};

namespace module {

	class Problem2 : public Module
	{
	public:
		Problem2();
		~Problem2();

		void OnUpdate(double deltaTime) override;
		void OnRender() override;
		void OnImguiRender() override;
		void CursorPosCallback(GLFWwindow* window, double xpos, double ypos) override;
		void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) override;
		void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) override;
		void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset) override;

		void GenerateMesh();
		void ClearTemperature(double initialTemperature);
		void TemperatureUpdate();
		void ColorUpdate();
		void PrintTemperature();
		void PrintMatrixK();
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

		std::vector<Vertex> m_Vertices;
		std::vector<Face> m_Faces;
		std::vector<unsigned int> m_Indices;

		ImGuiIO& m_IO;

		float m_scale;
		int n1;//网格径向数量
		int n2;//网格环向数量
		BoundaryType insideBoundary;
		BoundaryType outsideBoundary;
		double insideBoundaryValue;//固定的温度值
		double outsideBoundaryValue;//固定的流入热量值
	};
}
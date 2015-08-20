#pragma once

#ifndef GLUT_HELPER_H_
#define GLUT_HELPER_H_

//STL
#include <string>
#include <vector>
#include <iostream>

//GL
#include <GL\glew.h>
#include <gl\GL.h>
#include <GL\freeglut.h>

struct GLUTSettings
{
	int positionX;
	int positionY;
	int width;
	int height;
	std::string windowName;
	int GLVersionMajor;
	int GLVersionMinor;
};

struct FORMAT_RGBA_FLOAT
{
	float R;
	float G;
	float B;
	float A;
};

int I2_2_I2(const int row, const int col, const int width)
{
	return row * width + col;
}

void GLUT_HELPER_displayLoop(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	glutSwapBuffers();
}

int moving, beginx, beginy;
GLfloat anglex = 0;   /* in degrees */
GLfloat angley = 0;   /* in degrees */

void
motion(int x, int y)
{
	if (moving) {
		anglex = anglex + (x - beginx);
		angley = angley + (y - beginy);
		beginx = x;
		beginy = y;
		glutPostRedisplay();
	}

	glRotatef(angley, 1.0, 0.0, 0.0);
	glRotatef(-anglex, 0.0, 0.0, 1.0);
}

void mouseFunction(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		moving = 1;
		beginx = x;
		beginy = y;
	}
	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
		moving = 0;
	}
}

class GLUTHelper
{
public:

	GLUTHelper();

	void initWindow(int argc, char* argv[], const GLUTSettings& settings);

	void initCamera(float x, float y, float z, float lookAtX, float lookAtY, float lookAtZ);

	void setIdleFunc(void(*func)(void));

	void enterDisplayLoop(void(*func)(void));

	void renderTexture(const std::vector<FORMAT_RGBA_FLOAT>& values, int width, int height);

	void renderTexturePatch();

	void renderMesh();

private:

	bool m_isInitialised;
	bool m_textureInitialised;
	GLuint m_texture;
	GLuint m_textureQuad;
	GLuint m_textureQuadIndices;
};

//IMPLEMENTATION-------------------------------------------------------------------------------

GLUTHelper::GLUTHelper()
{
	m_isInitialised = false;
}

void GLUTHelper::setIdleFunc(void(*func)(void))
{
	glutIdleFunc(func);
}

void GLUTHelper::initCamera(float x, float y, float z, float lookAtX, float lookAtY, float lookAtZ)
{
	gluPerspective( /* field of view in degree */ 90.0,
		/* aspect ratio */ 1.0,
		/* Z near */ 1.0, /* Z far */ 500000.0);
	glMatrixMode(GL_MODELVIEW);

	gluLookAt(x, y, z,  /* eye is at (0,0,5) */
		lookAtX, lookAtY, lookAtZ,      /* center is at (0,0,0) */
		0.0, 1.0, 0.0);      /* up is in positive Y direction */
}

void GLUTHelper::enterDisplayLoop(void(*func)(void) = nullptr)
{
	if (func != nullptr)
	{
		glutDisplayFunc(func);
	}
	else
	{
		std::cout << "Entering Default Display Loop..." << std::endl;
		glutDisplayFunc(GLUT_HELPER_displayLoop);
	}

	glutMainLoop();
}

void GLUTHelper::initWindow(int argc, char* argv[], const GLUTSettings& settings)
{
	if (m_isInitialised)
	{
		std::cerr << "ERROR: GLUT & GLEW already initialised!" << std::endl;
		return;
	}

	//init GLUT
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_MULTISAMPLE);
	glutInitWindowPosition(settings.positionX, settings.positionY);
	glutInitWindowSize(settings.width, settings.height);
	glutCreateWindow(settings.windowName.c_str());
	glutInitContextVersion(settings.GLVersionMajor, settings.GLVersionMinor);
	glutMouseFunc(mouseFunction);
	glutMotionFunc(motion);

	//init GLEW
	GLenum status = glewInit();
	if (status != GLEW_OK)
	{
		std::cerr << "GLEW INIT ERROR: " << glewGetErrorString(status) << std::endl;
	}
	else
	{
		std::cout << "GLUT & GLEW initialised successfully. GLEW_VERSION: "
			<< glewGetString(GLEW_VERSION) << std::endl;
	}

	m_isInitialised = true;
}

void GLUTHelper::renderTexture(const std::vector<FORMAT_RGBA_FLOAT>& values, int width, int height)
{
	if (!m_textureInitialised)
	{
		//Coordinates for quad we wish to texture
		//e.g: http://stackoverflow.com/questions/18868844/draw-textured-quad-in-background-of-opengl-scene

		//Vertex & UV Coordinates
		GLfloat vertexUVCoords[] = { -1.0f, -1.0f, 0.0f,
			0.0f, 0.0f, 0.0f,
			1.0f, -1.0f, 0.0f,
			1.0f, 0.0f, 0.0f,
			1.0f, 1.0f, 0.0f,
			-1.0f, 1.0f, 0.0f,
			0.0f, 1.0f, 0.0f };

		GLubyte indices[] = { 0, 1, 2, 0, 2, 3 };

		//1. Generate VAO name
		glGenVertexArrays(1, &m_texture);
		//2. Bind VAO
		glBindVertexArray(m_texture);

		//3. Generate VBO
		glGenBuffers(1, &m_textureQuad);
		//4. Bind VBO
		glBindBuffer(GL_ARRAY_BUFFER, m_textureQuad);
		//5. Create & initialise VBO (using GL_STATIC_DRAW because we don't expect VBO to change)
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertexUVCoords), vertexUVCoords, GL_STATIC_DRAW);

		//6. Generate VBO Index buffer
		glGenBuffers(1, &m_textureQuadIndices);
		//7. Bind buffer
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_textureQuadIndices);
		//8. Create & initialise buffer
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);


	}
	else
	{
		//1. Update existing texture
	}
}

void GLUTHelper::renderTexturePatch()
{

}

void GLUTHelper::renderMesh()
{

}


#endif //GLUT_HELPER_H_
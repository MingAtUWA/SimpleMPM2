#include "PostProcessor_pcp.h"

#include "gif.h"

#include "GenerateAnimation.h"

GenerateAnimation::GenerateAnimation() :
	time_rcds(nullptr), time_rcd_num(0) {}
GenerateAnimation::~GenerateAnimation()
{
	if (time_rcds)
	{
		delete[] time_rcds;
		time_rcds = nullptr;
		time_rcd_num = 0;
	}
}

static void processInput(GLFWwindow *window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
}

static void framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
	GLsizei win_size;
	GLsizei win_padding;
	if (width < height)
	{
		win_size = width;
		win_padding = (height - width) / 2;
		glViewport(0, win_padding, win_size, win_size);
	}
	else
	{
		win_size = height;
		win_padding = (width - height) / 2;
		glViewport(win_padding, 0, win_size, win_size);
	}
}

static const size_t SCREEN_WIDTH = 600;
static const size_t SCREEN_HEIGHT = 600;

void reorder_buffer(unsigned char *RGBA_data, int width, int height)
{
	// RGBA data are 4 bytes long
	long *data = reinterpret_cast<long *>(RGBA_data);
	long *line1 = data;
	long *line2 = data + (height - 1) * width;
	long data_tmp;
	while (line1 < line2)
	{
		for (size_t i = 0; i < width; i++)
		{
			data_tmp = line1[i];
			line1[i] = line2[i];
			line2[i] = data_tmp;
		}
		line1 += width;
		line2 -= width;
	}
}

int GenerateAnimation::generate(double ani_time, double xl, double xu, double yl, double yu,
								const char *res_file_name, const char *gif_name)
{
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow *window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "model display", nullptr, nullptr);
	if (window == nullptr)
	{
		std::cout << "Failed to create GLFW window.\n";
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD.\n";
		glfwTerminate();
		return -2;
	}

	// init animation
	init(res_file_name);

	// animation_time
	if (time_rcd_num == 0)
		return -1; // no time record
	real_time = time_rcds[time_rcd_num-1].total_time - time_rcds[0].total_time;
	animation_time = ani_time;
	ani_real_ratio = animation_time / real_time;
	min_delay_real = 0.01 / ani_real_ratio;

	// gif output
	GifWriter gif_file;
	unsigned char *pixels_data;
	if (gif_name)
	{
		GifBegin(&gif_file, gif_name, SCREEN_WIDTH, SCREEN_HEIGHT, 1);
		pixels_data = new unsigned char[SCREEN_WIDTH * SCREEN_HEIGHT * 4];
	}

	bool render_new_frame = true; // control whether it is time to swap and render new frame
	bool reach_last_frame = false; // whether run out of frame to render
	size_t draw_frame_id = 0;
	double prev_time;
	cur_time_rcd_id = 0;
	while (!glfwWindowShouldClose(window))
	{
		processInput(window);

		if (render_new_frame)
		{
			prev_time = glfwGetTime();
			render_frame(xl, xu, yl, yu);
			render_new_frame = false;
			// output to gif file
			if (!reach_last_frame)
			{
				std::cout << "frame " << draw_frame_id++
						  << " real time " << time_rcds[cur_time_rcd_id].total_time
						  << " ani time "  << time_rcds[cur_time_rcd_id].total_substep_num << "\n";
				reach_last_frame = !find_next_frame();
				if (gif_name)
				{
					// read pixel from back buffer
					glReadPixels(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, pixels_data);
					reorder_buffer(pixels_data, SCREEN_WIDTH, SCREEN_HEIGHT);
					GifWriteFrame(&gif_file, pixels_data, SCREEN_WIDTH, SCREEN_HEIGHT, (delay_ani_100th ? delay_ani_100th : 1));
				}
			}
		}

		if (glfwGetTime() - prev_time >= delay_ani)
		{
			glfwSwapBuffers(window);
			render_new_frame = true;
		}

		glfwPollEvents();
	}

	if (gif_name)
	{
		GifEnd(&gif_file);
		delete[] pixels_data;
		pixels_data = nullptr;
	}

	glfwTerminate();
	return 0;
}


bool GenerateAnimation::find_next_frame(void)
{
	double prev_time = time_rcds[cur_time_rcd_id].total_time;
	size_t frame_id;
	double diff_time;
	for (frame_id = cur_time_rcd_id + 1; frame_id < time_rcd_num; ++frame_id)
	{
		diff_time = time_rcds[frame_id].total_time - prev_time;
		if (diff_time >= min_delay_real)
		{
			delay_ani = diff_time * ani_real_ratio;
			delay_ani_100th = unsigned short int(delay_ani * 100);
			break;
		}
	}
	
	if (frame_id >= time_rcd_num)
	{
		delay_ani = 0.25; // update every 0.25s
		delay_ani_100th = 25;
		return false;
	}

	cur_time_rcd_id = frame_id;
	return true;
}

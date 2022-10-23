#include "simple_interface.h"
#include <iostream>

int main(void)
{
	SimpleInterface detector;
	cv::Mat img1;
	cv::Mat img2;

	img1 = cv::imread("C:\\Users\\plonc\\OneDrive\\Documents\\src\\github_mio\\siftDetector\\uch010a.jpg");
	img2 = cv::imread("C:\\Users\\plonc\\OneDrive\\Documents\\src\\github_mio\\siftDetector\\uch010b.jpg");

	detector.addImg(img1, "case");
	detector.build();

	std::vector<Detection> detections;
	cv::Mat imdraw = detector.detect(img2, detections);

	std::cout << detections.size() << " detections" << std::endl;

	for (size_t i = 0; i < detections.size(); i++)
	{
		Detection& d = detections[i];
		std::cout << "  " << d.tx << " " << d.ty << " " << d.label << std::endl;
	}

	cv::imwrite("C:\\Users\\plonc\\OneDrive\\Documents\\src\\github_mio\\siftDetector\\imdraw.jpg", imdraw);

	return 0;
}
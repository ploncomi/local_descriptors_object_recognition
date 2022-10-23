#include "simple_interface.h"

SimpleInterface::SimpleInterface() : initialized_(false), trained_(false)
{

}

void SimpleInterface::addImg(cv::Mat img, std::string label_name)
{
	if (initialized_ == false)
	{
		rec_.fijaReconocedor(&sift_, &kdtree_, &hough_);
		//rec_.leeParametrosPropaga("C:\\Users\\plonc\\OneDrive\\Documents\\src\\github_mio\\siftDetector\\ParamSIFT.txt");
		initialized_ = true;
	}
	if (trained_ == true)
	{
		clear();
	}
	L_ImageGrayDouble imBN;
	imBN = img;
	L_ImageGrayUchar imGray;
	imGray = imBN;
	rec_.agregaImagenEntrenRobaObjetoDe(imBN);
	labels_.push_back(label_name);
}

void SimpleInterface::build()
{
	if (trained_ == true)
	{
		printf("SimpleInterface::build() Error: System training is done yet\n");
		return;
	}
	rec_.finalizaEntren();
	trained_ = true;
}

cv::Mat SimpleInterface::detect(cv::Mat img, std::vector<Detection>& detections)
{

	if (trained_ == false)
	{
		printf("SimpleInterface::detect() Error: System training is not done yet\n");
		return cv::Mat();
	}
	L_ImageGrayDouble imBN;
	imBN = img;
	rec_.procesaImagenesPrueba(imBN);

	//std::cout << rec_.calL.n << " calces" << std::endl;

	L_ShapeArray lins;
	L_ImageRGBUchar imuchar;
	rec_.dibujaTransformaciones(imuchar, 0, 1000, 1000, 1000);

	cv::Mat output;
	imuchar.copyTo(output);
	return output;
}

void SimpleInterface::clear()
{
	rec_.reseteaEntren();
	labels_.clear();
}


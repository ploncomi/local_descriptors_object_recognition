#ifndef SIMPLE_INTERFACE_H_
#define SIMPLE_INTERFACE_H_

#include "L_RecPuntInt/L_Descriptores.h"
#include "L_RecPuntInt/L_Fnes_Pato.h"
#include "L_RecPuntInt/L_GenPuntInt.h"
#include "L_RecPuntInt/L_RecPuntInt.h"

#include "L_RecPuntInt/siftImpl/L_CvSift.h"

#include <vector>
#include <string>

struct Detection
{
	double mxx; // Affine transform coefficient
	double mxy; // Affine transform coefficient
	double tx;  // Affine transform (translation in x)
	double myx; // Affine transform coefficient
	double myy; // Affine transform coefficient
	double ty;  // Affine transform (translation in y)
	double cx;  // Object's center
	double cy;  // Object's center
	double w;   // Width of bounding box
	double h;   // Height of bounding box
	std::string label;
};

class SimpleInterface
{
public:
	SimpleInterface();
	void addImg(cv::Mat img, std::string label_name);
	void build();
	cv::Mat detect(cv::Mat img, std::vector<Detection>& detections);
	void clear();

	bool initialized_;
	bool trained_;
	L_CvSift sift_;
	//L_KdTreeBBF kdtree_;
	L_Flann kdtree_;
	L_GenTransfAfinHough hough_;
	L_RecPuntInt rec_;

	std::vector<std::string> labels_;
};

#endif // SIMPLE_INTERFACE_H_


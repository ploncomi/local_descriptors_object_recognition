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
	double mxx;
	double mxy;
	double tx;
	double myx;
	double myy;
	double ty;
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
	L_KdTreeBBF kdtree_;
	L_GenTransfAfinHough hough_;
	L_RecPuntInt rec_;

	std::vector<std::string> labels_;
};

#endif // SIMPLE_INTERFACE_H_


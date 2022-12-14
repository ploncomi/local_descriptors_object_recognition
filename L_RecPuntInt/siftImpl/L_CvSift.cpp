#ifdef __COMPAT_CVMAT__
#include "L_CvSift.h"

// Works with opencv 4.4+

bool L_CvSift::calcPuntInt_act()
{
	return false;
}

bool L_CvSift::calcDirecc_act()
{
	return false;
}

bool L_CvSift::calcDescr_act()
{
	return false;
}

bool L_CvSift::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	L_ImageGrayUchar imBN;
	imBN = im;
	cv::Mat mat;
	imBN.copyTo(mat);

	// Define features detector
	cv::Ptr<cv::SiftFeatureDetector> detector = cv::SiftFeatureDetector::create(150, 6);
	// Detect the keypoints
	std::vector<cv::KeyPoint> keypoints;

	detector->detect(mat, keypoints);

	cv::Mat descriptors;
	cv::Ptr<cv::SiftDescriptorExtractor> extractor = cv::SiftDescriptorExtractor::create(150, 6);
	extractor->compute(mat, keypoints, descriptors);

	L_Descriptor d;
	d.tipoP = tipoP;
	d.tipoD = tipoD;
	d.nobj = nobj;
	d.nt = 128;
	d.imorig = &im;
	d.imorig_lx = im.lx;
	d.imorig_ly = im.ly;
	d.hasDepth = 0;
	d.refGen = this;

	if (int(keypoints.size()) != descriptors.rows)
	{
		printf("Error in L_CvSift::calcDescrPuntInt(): Wrong keypoints / descriptors size\n");
		return false;
	}

	if (descriptors.rows>0 && descriptors.cols != d.nt)
	{
		printf("Error in L_CvSift::calcDescrPuntInt(): Wrong descriptor length\n");
		return false;
	}

	for (size_t i = 0; i < keypoints.size(); i++)
	{
		d.x0 = keypoints[i].pt.x;
		d.y0 = keypoints[i].pt.y;
		d.sigma0 = keypoints[i].size;
		d.ang = -keypoints[i].angle * M_PI / 180;
		d.intens = fabs(keypoints[i].response);  // Aca no estoy considerando el signo ->laplacian
		d.signo = keypoints[i].response > 0 ? 1 : -1;
		d.id = ultId++;
		if (descriptors.type() == CV_8UC1)
		{
			for (int j = 0; j < descriptors.cols; j++)
				d.vector[j] = descriptors.at<unsigned char>(i, j) / 255.0;
		}
		else if (descriptors.type() == CV_32FC1)
		{
			for (int j = 0; j < descriptors.cols; j++)
				d.vector[j] = descriptors.at<float>(i, j) / 255.0;
		}
		else
			printf("Bad description format");

		desL.push_back(d);
	}

	return true;
}

#endif // __COMPAT_CVMAT__

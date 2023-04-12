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

	for (L_TransfPuntInt2DNodoPtr* ptr = rec_.traL.root; ptr != NULL; ptr = ptr->sig)
	{
		Detection det;
		L_TransfAfinPI2D* tr = (L_TransfAfinPI2D*)ptr->c;
		det.mxx = tr->m11;
		det.mxy = tr->m12;
		det.tx = tr->tx;
		det.myx = tr->m21;
		det.myy = tr->m22;
		det.ty = tr->ty;
		double lx = ptr->c->calL.root->c.dRef->imorig_lx;
		double ly = ptr->c->calL.root->c.dRef->imorig_ly;
		det.cx = det.mxx * lx/2 + det.mxy * ly / 2 + det.tx;
		det.cy = det.myx * lx / 2 + det.myy * ly / 2 + det.ty;
		double p1x = det.mxx * 0 + det.mxy * 0 + det.tx;
		double p1y = det.myx * 0 + det.myy * 0 + det.ty;
		double p2x = det.mxx * lx + det.mxy * 0 + det.tx;
		double p2y = det.myx * lx + det.myy * 0 + det.ty;
		double p3x = det.mxx * lx + det.mxy * ly + det.tx;
		double p3y = det.myx * lx + det.myy * ly + det.ty;
		double p4x = det.mxx * 0 + det.mxy * ly + det.tx;
		double p4y = det.myx * 0 + det.myy * ly + det.ty;
		double x1 = std::min(std::min(p1x, p2x), std::min(p3x, p4x));
		double y1 = std::min(std::min(p1y, p2y), std::min(p3y, p4y));
		double x2 = std::max(std::max(p1x, p2x), std::max(p3x, p4x));
		double y2 = std::max(std::max(p1y, p2y), std::max(p3y, p4y));
		det.w = x2 - x1;
		det.h = y2 - y1;
		det.label = labels_[ptr->c->calL.root->c.dRef->nobj];
		detections.push_back(det);
	}

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


#pragma once

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <random>
#include <time.h>
#include <fstream>

using namespace std;
using namespace cv;

struct ellipse_data
{
    double a, b, angle;
    double x0, y0;
    
    ellipse_data()
    {
        a = 0; b = 0; angle = 0; x0 = 0; y0 = 0;
    }
    
    ellipse_data(double _x0, double _y0, double _a, double _b, double th)
    {
        x0 = _x0; y0 = _y0; a = _a; b = _b; angle = th;
    }
    
};

struct sort_by_score
{
    bool operator()(std::pair<ellipse_data, double>const& a, std::pair<ellipse_data, double>const& b)
    {
        return a.second > b.second;
    }
};

int detect_ellipse(vector<pair<ellipse_data, double>>& result,const vector<Point_<int>>& edges, int n_best);

vector<vector<int>> distances_array(const vector<Point_<int>>& points);

vector<Point_<int>> find_close_points(const vector<vector<int>>& points, int min_value, int max_value);

vector<double> get_angular_array(const vector<Point_<int>>& edges, const vector<Point_<int>>& close_points);

void shuffle_vector(vector<Point_<int>>& points, int max_N);

vector<double> accumarray(const vector<int>& indexes, const vector<double>& weights,  int sz_1);

vector<Point_<int>> find_edges(Mat& src);

vector<double> get_weights_matrix(const Mat& src, const vector<Point_<int>>& points);

vector<double> convolve(vector<double>& accum, const vector<double>& kernel, int zero_lengths);

pair<int, double> get_best_score_and_idx(const vector<double>& accum);

vector<double> get_kernel(int extension, double str_dev);

vector<Point_<int>> choose_right_angles(const vector<Point_<int>>& points, const vector<double>& angles,
                                        double min_a, double max_a, bool mode_high);

int get_source_img(Mat& src,string& input_img, string& output_img, int argc, char** argv);

vector<Point_<int>> transform_cast(const vector<Point2d>& array);

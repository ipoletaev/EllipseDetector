#include "detector.h"

/*
 -- Copyright 2015-2016  Author: Poletaev Igor --
 
 Input arguments:
 --------
 inp_img
 - One-channel input image (greyscale or binary)
 out_img
 - name of output image with labeled ellipse
 
 Parameters:
 --------
 - minMajorAxis: Minimal length of major axis accepted.
 - maxMajorAxis: Maximal length of major axis accepted.
 - rotation, rotationSpan: Specification of restriction on the angle of the major axis in degrees.
 If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,
 rotation+rotationSpan] are accepted.
 - minAspectRatio: Minimal aspect ratio of an ellipse (in (0,1))
 - randomize: Subsampling of all possible point pairs. Instead of examining all N*N pairs, runs
 only on N*randomize pairs. If 0, randomization is turned off.
 - numBest: Top numBest to return
 - uniformWeights: Used to prefer some points over others. If false, accumulator points are weighted
 by their grey intensity in the image. If true, the input image is regarded as binary.
 - smoothStddev: In order to provide more stability of the solution, the accumulator is convolved with
 a gaussian kernel. This parameter specifies its standard deviation in pixels.
 
 Return values:
 --------
 Returns a matrix of best fits. Each row (there are params.numBest of them) contains six elements:
 [[x0 y0 a b alpha] score]] being the center of the ellipse, its major and minor axis, its angle in degrees and score.
 */


// **** parameters ***** //
double min_major_axis = 160;
double max_major_axis = 240;
double rotation = 0;
double rotation_span = 0;
double min_aspect_ratio = 0.1;
double smooth_std_dev = 1.0;
double eps = 0.0001;
bool uniform_weights = true;
int randomize = 1;
// ********************* //

int detect_ellipse(vector<pair<ellipse_data, double>>& result,const vector<Point_<int>>& edges, int n_best)
{
    // ******** intializations ********* //
    if (rotation_span > 90) rotation_span = 90;
    vector<double> gaussian_kernel = get_kernel(6, smooth_std_dev);
    size_t N = edges.size();
    cout << "Possible major axis: " << N * N << endl;
    
    vector<vector<int>> distances = distances_array(edges);
    vector<Point_<int>> close_points = find_close_points(distances, min_major_axis, max_major_axis);
    N = close_points.size();
    cout << "Possible major axis after distances costraints: " << N << endl;
    // ********************************* //
    
    // ********* constraints *********** //
    if (rotation_span > 0)
    {
        //angle constrains
        double low_tan = tan(rotation - rotation_span);
        double high_tan = tan(rotation + rotation_span);
        vector<double> tangents = get_angular_array(edges, close_points);
        close_points = choose_right_angles(close_points, tangents, low_tan, high_tan, (high_tan > low_tan));
        N = close_points.size();
        cout << "Possible major axis after angular costraints: " << N << endl;
    }
    else
    {
        cout << "Angular constraints aren't using ... " << endl;
    }
    if (randomize > 0)
    {
        shuffle_vector(close_points, (int)min(edges.size() * randomize, N));
        cout << "Get " << close_points.size() << " random points to process ..." << endl;
    }
    // ********************************* //
    
    // ********** main flow ************ //
    int x1,y1,x2,y2, current_dist, minor_propositions;
    double x0, y0, a_sq, current_cos_tau;
    for (const Point_<int>& p: close_points)
    {
        x1 = edges[p.x].x; x2 = edges[p.y].x;
        y1 = edges[p.x].y; y2 = edges[p.y].y;
        
        //compute center & major axis
        x0 = (x1 + x2) / 2.0;
        y0 = (y1 + y2) / 2.0;
        a_sq = distances[p.x][p.y] / 4.0;
        
        vector<Point_<int>> current_points;
        vector<int> sq_distances;
        vector<double> tau_cosine;
        for (int i = 0; i < edges.size(); i++)
        {
            current_dist = pow(edges[i].x - x0, 2) + pow(edges[i].y - y0, 2);
            if (current_dist <= a_sq)
            {
                minor_propositions = pow(edges[i].x - x2, 2) + pow(edges[i].y - y2, 2);
                current_cos_tau = (a_sq + current_dist - minor_propositions)/(2.0 * sqrt(a_sq * current_dist));
                
                if (current_cos_tau > 1.0) current_cos_tau = 1.0;
                else if (current_cos_tau < -1.0) current_cos_tau = -1.0;
                
                tau_cosine.push_back(current_cos_tau);
                sq_distances.push_back(current_dist); // otherwise the formulae in paper does not work
                current_points.push_back(Point_<int>(edges[i].x, edges[i].y));
            }
        }
        
        vector<int> indexes(tau_cosine.size(), 0);
        for (int i = 0; i < tau_cosine.size(); i++)
        {
            indexes[i] = (int)(sqrt(a_sq * (1 - pow(tau_cosine[i], 2)) * sq_distances[i] /
                                    (a_sq - sq_distances[i] * pow(tau_cosine[i], 2) + eps)));
        }
        vector<double> weights(current_points.size(), 1.0); //= get_weights_matrix(src, current_points);
        
        // voiting
        vector<double> accum = accumarray(indexes, weights, max_major_axis);
        
        //a bit of smoothing and finding the most busy bin
        accum = convolve(accum, gaussian_kernel, (int)(sqrt(a_sq) * min_aspect_ratio));
        pair<int, double> score = get_best_score_and_idx(accum);
        int idx = score.first;
        double best_score = score.second;
        result.push_back(make_pair(ellipse_data(x0, y0, sqrt(a_sq), idx, atan((y1 - y2) / (double)(x1 - x2))), best_score));
    }
    // ********************************* //
    
    sort(result.begin(), result.end(), sort_by_score());
    result.resize(n_best);
    
    return 1; // success
}

vector<vector<int>> distances_array(const vector<Point_<int>>& points)
{
    size_t max_idx = points.size();
    vector<vector<int>> answer(max_idx, vector<int>(max_idx, 0.0));
    for (int i = 0; i < max_idx; i ++)
    {
        for (int j = 0; j < i; j ++)
        {
            answer[i][j] = answer[j][i] = pow(points[i].x - points[j].x, 2) + pow(points[i].y - points[j].y, 2);
        }
    }
    
    return answer;
}

vector<Point_<int>> find_close_points(const vector<vector<int>>& points, int min_value, int max_value)
{
    size_t max_idx = points.size();
    vector<Point_<int>> answer;
    for (int i = 0; i < max_idx; i++)
    {
        for (int j = 0; j < max_idx; j++)
        {
            long int mult = points[i][j];
            if (mult >= min_value * min_value && mult <= max_value * max_value && i < j)
            {
                answer.push_back(Point_<int>(i, j));
            }
        }
    }
    
    return answer;
}

vector<double> get_angular_array(const vector<Point_<int>>& edges, const vector<Point_<int>>& close_points)
{
    vector<double> answer(close_points.size(), 0.0);
    int k = 0;
    
    for (const Point_<int>& point: close_points)
    {
        answer[k++] = (edges[point.x].y - edges[point.y].y) / (double)(edges[point.x].x - edges[point.y].x); //tan(a)
    }
    
    return answer;
}

vector<Point_<int>> choose_right_angles(const vector<Point_<int>>& points, const vector<double>& angles,
                                      double min_a, double max_a, bool mode_high)
{
    vector<Point_<int>> clear_close_points;
    for (int i = 0; i < angles.size(); i++)
    {
        if (mode_high)
        {
            if (angles[i] > min_a && angles[i] < max_a) clear_close_points.push_back(points[i]);
        }
        else
        {
            if (angles[i] > min_a || angles[i] < max_a) clear_close_points.push_back(points[i]);
        }
    }
    
    return clear_close_points;
}

void shuffle_vector(vector<Point_<int>>& points, int max_N)
{
    if (max_N > 1)
    {
        size_t i;
        //srand(time(NULL));
        for (i = 0; i < max_N - 1; i++)
        {
            size_t j = i + rand() / (RAND_MAX / (max_N - i) + 1);
            swap(points[i], points[j]);
        }
    }
    
    points.resize(max_N);
}

vector<double> accumarray(const vector<int>& indexes, const vector<double>& weights,  int sz_1)
{
    
    /*
     -Find out how many unique indices there are in subs. Each unique index defines a bin in the output array. The maximum index value in subs determines the size of the output array.
     -Find out how many times each index is repeated.
     -This determines how many elements of vals are going to be accumulated at each bin in the output array.
     -Create an output array. The output array is of size max(subs) or of size sz.
     -Accumulate the entries in vals into bins using the values of the indices in subs and apply fun to the entries in each bin.
     -Fill the values in the output for positions not referred to by subs. Default fill value is zero; use fillval to set a different value.
     */
    
    if (indexes.size() == weights.size())
    {
        vector<double> answer(sz_1, 0.0);
        for (int i = 0; i < sz_1; i++)
        {
            vector<int> temp_idx;
            double sum = 0.0;
            
            for (int j = 0; j < indexes.size(); j++)
            {
                if (indexes[j] == (i + 1)) temp_idx.push_back(j);
            }
            if (!temp_idx.empty())
            {
                for (int j = 0; j < temp_idx.size(); j++)
                {
                    sum += weights[temp_idx[j]];
                }
                
                answer[i] = sum;
            }
            else
            {
                answer[i] = 0;
            }
        }
        
        return answer;
    }
    else
    {
        cerr << "Dimensions of the input vectors are mismatched!" << endl;
        return vector<double>(0);
    }
}

vector<Point_<int>> find_edges(Mat& src)
{
    //Canny(src, src, 50, 250, 3);
    //imwrite("/Users/Igor/Desktop/test.png", src);
    vector<Point_<int>> answer;
    for (int i = 0; i < src.rows; i ++)
    {
        for (int j = 0; j < src.cols; j ++)
        {
            if ((int)src.at<uchar>(i, j) > 0)
            {
                answer.push_back(Point_<int>(i, j));
            }
        }
    }
    
    return answer;
}

vector<double> get_weights_matrix(const Mat& src, const vector<Point_<int>>& points)
{
    
    if (uniform_weights)
    {
        return vector<double>(points.size(), 1.0);
    }
    
    vector<double> weights(points.size(), 0.0);
    for (int i = 0; i < points.size(); i++)
    {
        if ((int)src.at<uchar>(points[i].x, points[i].y) > 0) weights[i] = (int)src.at<uchar>(points[i].x, points[i].y) / 255.0; // we suggest to use grayscale images only
    }
    
    return weights;
}

vector<double> convolve(vector<double>& accum, const vector<double>& kernel, int zero_lengths)
{
    int const singnal_size = accum.size();
    int const kernel_size= kernel.size();
    int const n  = singnal_size + kernel_size - 1;
    vector<double> full_answer(n, 0.0);
    
    for(auto i(0); i < n; ++i)
    {
        int const j_min = (i >= kernel_size - 1)? i - (kernel_size - 1) : 0;
        int const j_max = (i < singnal_size - 1)? i : singnal_size - 1;
        for(auto j = j_min; j <= j_max; ++j)
        {
            full_answer[i] += (accum[j] * kernel[i - j]);
        }
    }
    
    int shift = floor(kernel_size / 2.0);
    vector<double> answer(full_answer.begin() + shift, full_answer.begin() + singnal_size + shift);
    for (int i = 0; i < zero_lengths; i++)
    {
        answer[i] = 0.0;
    }
    
    return answer;
}

pair<int, double> get_best_score_and_idx(const vector<double>& accum)
{
    double max = accum[0];
    int idx = 0;
    for (int i = 0; i < accum.size(); i++)
    {
        if (accum[i] > max)
        {
            max = accum[i]; idx = i;
        }
    }
    
    return make_pair(idx, max);
}

vector<double> get_kernel(int extension, double str_dev)
{
    
    // get lowpass Gaussian 1D filter of size = <extension> with standard deviation = <std_dev>
    Mat q = getGaussianKernel(extension, str_dev);
    vector<double> answer(extension, 0.0);
    for (int a = 0; a < extension; a++)
    {
        answer[a] = q.at<double>(a);
    }
    return answer;
}

int get_source_img(Mat& src,string& input_img, string& output_img, int argc, char** argv)
{
    if (argc < 2)
    {
        cout << "Usage: <input file name> [<output file name> (default: \"ellipses.jpg\")"
        << "[<maximum number of ellipces to find (default: 1)>]]" << endl;
        return -1;
    }
    
    string inp_img = argv[1];
    if (argc > 2)
    {
        output_img = argv[2];
    }
    else
    {
        output_img = "out.png";
    }
    
    src = imread(inp_img, 0);
    if(! src.data )
    {
        cout <<  "Could not open or find the image" << std::endl;
        return -1;
    }
    
    input_img = inp_img;
    return 1;
}

vector<Point_<int>> transform_cast(const vector<Point2d>& array)
{
    vector<Point_<int>> answer(array.size());
    for (int i = 0; i < array.size(); i++)
    {
        answer[i].x = floor(array[i].x);
        answer[i].y = floor(array[i].y);
    }
    return answer;
}


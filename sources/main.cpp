#include "detector.h"
#include <iostream>
#include <vector>
#include <iterator>
#include <fstream>
#include <algorithm>

using namespace std;

int n_best = 5;

int main(int argc, char** argv)
{
    srand(time(NULL));
    Mat source;
    string inp_img, out_img;
    
    if (get_source_img(source, inp_img, out_img, argc, argv) == 1)
    {
        Mat ellipses_draw = imread(argv[1], 1);
        vector<pair<ellipse_data, double>> best_fit(3);
        
        if (detect_ellipse(best_fit, find_edges(source), n_best))
        {
            for (int k = 0; k < n_best; k ++)
            {
                ellipse_data elp = best_fit[k].first;
                if ((elp.x0 >= 0.0) && (elp.y0 >= 0.0))
                {
                    cout << "Found: (x_0 = " << elp.x0 << ", y_0 = " << elp.y0 << ", a = " << elp.a
                    << ", b = " << elp.b << ", angle = " << elp.angle << ") with score = " << best_fit[k].second << endl;
                
                    // draw an ellipse via means of openCV
                    ellipse(ellipses_draw,
                        Point_<int>((int)elp.y0, (int)elp.x0),
                        Size_<double>(elp.a, elp.b),
                        elp.angle * 180 / CV_PI,
                        0, 360, Scalar(0, 0, 255), 1, LINE_AA);
                }
            }
            
        imwrite(out_img, ellipses_draw);
            
        }
    }

    return 0;
    
}

// usage example:
ellipse_data fit_ellipse(const vector<Point2d>& data)
{
    vector<pair<ellipse_data, double>> best_fit(1);
    if (detect_ellipse(best_fit, transform_cast(data), 1))
    {
        ellipse_data elp = best_fit[0].first;
        if ((elp.x0 >= 0.0) && (elp.y0 >= 0.0) && best_fit[0].second > 0)
        {
            return elp;
        }
    }
    return ellipse_data();
}

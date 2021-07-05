#include "ImageDisplaceLib.h"
using namespace std;

void testing(){
    vector<double> star = {54.12015333986315, 12.65005233482232, 0.25993125876773404, 0.8111479866192958};
    vector<double> ray = {0, 0};
    pair<double, double> res;

    res = get_displace_of_ray_by_single_star(star, ray);

    cout << res.first << " " << metres_to_mks(res.second);


}


void get_path(int x, int y, string code){
    bool release = true;
    if (!release) {
        testing();

    }

    std::uniform_real_distribution<> urd(0, 10);


    vector<vector<double>> stars = get_stars_set("C:\\RW\\stars.txt", 3597437, false);
    vector<vector<double>> velocities (stars.size());
    vector<vector<double>> corona_velocities (stars.size());

    corona_velocities = get_random_3d_velocity(stars.size(), 200, 250);
    velocities = get_random_3d_velocity(stars.size(), 15, 25);
    for (long i = 0; i < stars.size(); i++) {
        if (stars[i][4] > 0) velocities[i][0] += 225000;
        if (stars[i][4] == 0) velocities[i] = corona_velocities[i];
    }

    pair<int, int> time_config = {10000, 2000};
    long long full_duration = sec_in_year * time_config.first;
    long long interval = sec_in_day * time_config.second;
    double interval_double = interval;

    string path_to_file;
    path_to_file = "C:\\RW\\path" + code + ".txt";
    ofstream fout(path_to_file);




    pair<double, double> image;
    vector<double> ray = {x * pc, y * pc};
    fout << "Path for ray(in pc) " << ray[0] / pc << " " << ray[1] / pc << endl;
    fout << "Points calculated each " << time_config.second << " days for " << time_config.first << " years." << endl;

    for (long long time_passed = 0; time_passed < full_duration; time_passed += interval){
        cout << "Calculating day " << time_passed / interval * time_config.second << endl;
        image = get_full_displace(stars, ray);
        fout << metres_to_mks(image.first) << " " << metres_to_mks(image.second) << endl;

        for (long i = 0; i < stars.size(); i++) for (int j = 0; j < 3; j++) {
                stars[i][j] += (velocities[i][j] * interval_double) / pc;
            }
    }


    fout.close();
}

int main(){
    //get_path(1, 1, "long");
    double dist_btw_rays = tan(0.2 * 1000000 / mks) * full_H;

    //get_path(5, 11, "_pair11");
    //get_path(5, 11 - dist_btw_rays, "_pair12");
    get_path(5 + dist_btw_rays, 11, "_pair13");

    return 0;


    int n = 7;
    for (int i = -n; i <= n; i++){
        for (int j = -n; j <= n; j++){
            get_path(i, j, to_string(i) + "_" + to_string(j));
        }
    }
}


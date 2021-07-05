#include "ImageDisplaceLib.h"

double mks_to_metres(double angle){
    double metres;
    metres = (angle * full_H * pc) / mks;
    return metres;
}

double metres_to_mks(double metres){
    double angle;
    angle = (metres * mks) / (full_H * pc);
    return angle;
}

void print(std::vector<double>& lst, bool new_line, bool error_stream) {
    for (double x : lst) {
        if (error_stream == 0) std::cout << x << " ";
        if (error_stream == 1) std::cerr << x << " ";
    }
    if (new_line) {
        if (error_stream == 0) std::cout << std::endl;
        if (error_stream == 1) std::cerr << std::endl;
    }

}


std::pair<double, double> get_displace_of_ray_by_single_star(std::vector<double> &star, std::vector<double> &ray) {
    // In vector star distances given in pc, mass in m0. Ray coord given in metres

    double temp_x = ray[0], temp_y = ray[1];
    double distance, alpha;
    std::pair<double, double> displace;

    temp_x *= (star[2] / full_H);
    temp_y *= (star[2] / full_H);
    temp_x -= star[0] * pc;
    temp_y -= star[1] * pc;

    distance = sqrt(temp_x * temp_x + temp_y * temp_y);

    if (distance < 1) {
        std::cerr << "WARNING: distance was less then 1 meter, it has been set to 100pc." << std::endl;
        std::cerr << "Ray and star were ";
        print(ray, false, true);
        print(star, true, true);

        distance = 100 * pc;
    }
    alpha = (4 * G * star[3] * m0) / (c * c * distance);

    double multiplier;
    multiplier = (alpha * full_H * pc) / distance;

    displace.first = temp_x * multiplier;
    displace.second = temp_y * multiplier;

    return displace;
}

std::pair<double, double> get_full_displace(std::vector<std::vector<double>> &stars, std::vector<double> &ray){
    std::pair<double, double> temp, displace = {0, 0};

    for (auto & star : stars){
        temp = get_displace_of_ray_by_single_star(star, ray);
        displace.first += temp.first;
        displace.second += temp.second;
    }

    return displace;
}

double get_correlation_of_values(std::vector<double> &v1, std::vector<double> &v2){
    if (v1.size() != v2.size()) return 0;

    auto n = v1.size();
    double mean_value_v1 = 0, mean_value_v2 = 0;
    double covariance = 0, deviation1 = 0, deviation2 = 0;

    for (int i = 0; i < n; i++){
        mean_value_v1 += v1[i];
        mean_value_v2 += v2[i];
    }

    mean_value_v1 /= n;
    mean_value_v2 /= n;

    for (int i = 0; i < n; i++){
        covariance += (v1[i] - mean_value_v1) * (v2[i] - mean_value_v2);
        deviation1 += (v1[i] - mean_value_v1) * (v1[i] - mean_value_v1);
        deviation2 += (v2[i] - mean_value_v2) * (v2[i] - mean_value_v2);
    }

    if (deviation1 * deviation2 == 0) return 1;

    return (covariance / sqrt(deviation1 * deviation2));
}

std::vector<std::vector<double>> get_random_3d_velocity(long amount, double v_min, double v_max){
    std::vector<std::vector<double>> velocity(amount);
    std::vector<double> temp(3, 0);


    std::mt19937_64 gen(time(nullptr));
    std::uniform_real_distribution<> urd(-1, 1);
    std::uniform_real_distribution<> urd_max(v_min, v_max);
    double Vsum, V;

    for (long i = 0; i < amount; i++){
        temp[0] = urd(gen);
        temp[1] = urd(gen);
        temp[2] = urd(gen);
        Vsum = 1000 * urd_max(gen);
        V = sqrt(temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2]);
        temp[0] = temp[0] * Vsum / V;
        temp[1] = temp[1] * Vsum / V;
        temp[2] = temp[2] * Vsum / V;
        velocity[i] = temp;
    }

    return velocity;

}

std::vector<std::vector<double>> get_stars_set(std::string path, long n, bool generate_new){
    std::cout << "Getting stars' data..." << std::endl;
    if (generate_new){
        // Нужно разобраться с тонкостями работы system(command) и оптимизировать питоний код, может перенести сюда
        // По итогу этого ифа нужно оставлять в path путь на файл с звездами, в идеале еще создавать переменную
        // с их количеством, но лучше это делать не как на момент написания коммента - через отдельный файлик
        // Можно попробовать сразу звезды передавать, но что-то я сомневаюсь в простоте реализации такого варианта
        // system("python C:\\RW\\star_gen.py"); - эта штука работает, но вроде как адрес должен быть const
    }
    std::vector<double> single_star(5, 1); // последний элемент - тип - 1=disk, 2=flat, 0=corona
    char population_type;
    std::vector<std::vector<double>> stars(n);
    std::ifstream fin(path);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < 4; j++) fin >> single_star[j];
        fin >> population_type;
        if (population_type == 'd') single_star[4] = 1;
        if (population_type == 'f') single_star[4] = 2;
        if (population_type == 'c') single_star[4] = 0;
        stars[i] = single_star;
    }

    fin.close();

    return stars;
}
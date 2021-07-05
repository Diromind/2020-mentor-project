#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <random>
#include <ctime>
#include <fstream>

#ifndef CPP_PROJECTS_IMAGEDISPLACELIB_H
#define CPP_PROJECTS_IMAGEDISPLACELIB_H

const double pc = 30856775810000000.0; // метров
const double full_H = 390000.0; // парсеки
const double G = 6.674 / 1e+11; // в СИ
const double c = 299792458.0; // тоже в СИ
const double m0 = 1.989 * 1e+30; // масса солнца в кг
const double mks = 206265806000.0;  // в одном радиане микросекунд
const long long sec_in_day = 86400;
const long long sec_in_year = 31536000;



std::pair<double, double> get_displace_of_ray_by_single_star(std::vector<double> &star, std::vector<double> &ray);

std::pair<double, double> get_full_displace(std::vector<std::vector<double>> &stars, std::vector<double> &ray);



double get_correlation_of_values(std::vector<double> &v1, std::vector<double> &v2);



std::vector<std::vector<double>> get_stars_set(std::string path, long n, bool generate_new);

std::vector<std::vector<double>> get_random_3d_velocity(long amount, double v_min, double v_max);



void print(std::vector<double>& lst, bool new_line = true, bool error_stream = false);

double mks_to_metres(double angle);

double metres_to_mks(double metres);

#endif

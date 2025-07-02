#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <random>
#include <chrono>
#include <numeric> 

void RandomfillMatrix(std::vector<std::vector<int>>& matrix) {
    int N = matrix.size();
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            matrix[i][j] = rand() % (50 - 1 + 1) + 2;  
            matrix[j][i] = matrix[i][j];                 
        }
        matrix[i][i] = 0; 
    }
}

void FillMatrix(std::vector<std::vector<int>>& matrix) {
    int N = matrix.size();
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = i + 1; j < N; ++j)
        {
            int r;
            std::cin >> r;
            if (r < 0)
            {
                RandomfillMatrix(matrix);
                return;

            }
            matrix[i][j] = r;
            matrix[j][i] = matrix[i][j]; 
        }
    }
}
void PrintMatrix(std::vector<std::vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "___________________________"<<"\n";
}
void PrintVec(std::vector<int> vec){
    for (const auto& item : vec) {
        std::cout << item << " ";
    }
}
int LenCount(std::vector<int> route, std::vector<std::vector<double>>& matrix){
    int first_city = route[0];
    int start_city = first_city;
    int total = 0;
    for (auto it = route.begin() + 1; it != route.end(); ++it) {
        total+= matrix[start_city][*it];
        start_city = *it;
    }
    return total+ matrix[start_city][first_city];
}
void twoOptSwap(std::vector<int>& path, int i, int k) {
    while (i < k) {
        std::swap(path[i], path[k]);
        i++;
        k--;
    }
}
std::vector<int> createVector(int n) {
    std::vector<int> vec(n);

    // Заполняем числами от 0 до n-1
    std::iota(vec.begin(), vec.end(), 0);

    // Перемешиваем
    std::shuffle(vec.begin(), vec.end(), std::mt19937{ std::random_device{}() });

    return vec;
}

std::vector<int> hillClimbingTSP(int max_iterations, std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<int> currentPath = createVector(n);
    int currentDistance = LenCount(currentPath, matrix);

    for (int iter = 0; iter < max_iterations; iter++) {
        bool improved = false;
        std::vector<int> bestNeighbor = currentPath;
        int bestNeighborDistance = currentDistance;
        // Перебираем всех соседей (Swap двух городов)
        for (int i = 1; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                std::vector<int> neighbor = currentPath;
                std::swap(neighbor[i], neighbor[j]); // Меняем два города
                int neighborDistance = LenCount(neighbor, matrix);

                // Если нашли лучшее решение, переходим на него
                if (neighborDistance < bestNeighborDistance) {
                    bestNeighbor = neighbor;
                    bestNeighborDistance = neighborDistance;
                    improved = true;
                    currentPath = bestNeighbor;
                    currentDistance = bestNeighborDistance;
                    //PrintVec(bestNeighbor);
                    //std::cout << std::endl << "score " << bestNeighborDistance << std::endl;
                }
            }
        }
        break;
        // Проверяем всех соседей через 2-opt
        /*for (int i = 1; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                std::vector<int> neighbor = currentPath;
                twoOptSwap(neighbor, i, j);
                int neighborDistance = LenCount(neighbor, matrix);

                if (neighborDistance < bestNeighborDistance) {
                    bestNeighbor = neighbor;
                    bestNeighborDistance = neighborDistance;
                    improved = true;
                }
            }
        }*/

        // Если нашли улучшение, переходим к лучшему соседу
        
    }


    return currentPath;
}
std::vector<int> hillstohasticClimbingTSP(int max_iterations, std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<int> currentPath = createVector(n);
    int currentDistance = LenCount(currentPath, matrix);

    // Инициализация генератора случайных чисел
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    const double acceptance_probability = 0.1; // 10% шанс принять худшее решение

    for (int iter = 0; iter < max_iterations; iter++) {
        bool improved = false;
        std::vector<int> bestNeighbor = currentPath;
        int bestNeighborDistance = currentDistance;

        // Перебираем всех соседей (Swap двух городов)
        for (int i = 1; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                std::vector<int> neighbor = currentPath;
                std::swap(neighbor[i], neighbor[j]);
                int neighborDistance = LenCount(neighbor, matrix);

                // Если нашли лучшее решение или случайно принимаем худшее
                if (neighborDistance < bestNeighborDistance ||
                    (distribution(generator) < acceptance_probability)) {

                    bestNeighbor = neighbor;
                    bestNeighborDistance = neighborDistance;
                    improved = true;
                    currentPath = bestNeighbor;
                    currentDistance = bestNeighborDistance;
                    //PrintVec(bestNeighbor);
                    //std::cout << std::endl << "score " << bestNeighborDistance << std::endl;
                }
            }
        }

        if (!improved) {
            break; // Локальный оптимум достигнут
        }
    }

    return currentPath;
}
std::vector<int> hillsteepestClimbingTSP(int max_iterations, std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<int> currentPath = createVector(n);
    int currentDistance = LenCount(currentPath,matrix);

    for (int iter = 0; iter < max_iterations; iter++) {
        bool improved = false;
        std::vector<int> bestNeighbor = currentPath;
        int bestNeighborDistance = currentDistance;
        // Перебираем всех соседей (Swap двух городов)
        for (int i = 1; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                std::vector<int> neighbor = currentPath;
                std::swap(neighbor[i], neighbor[j]); // Меняем два города
                int neighborDistance = LenCount(neighbor,matrix);

                // Если нашли лучшее решение, переходим на него
                if (neighborDistance < bestNeighborDistance) {
                    bestNeighbor = neighbor;
                    bestNeighborDistance = neighborDistance;
                    improved = true;
                }
            }
        }
        // Проверяем всех соседей через 2-opt
        /*for (int i = 1; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                std::vector<int> neighbor = currentPath;
                twoOptSwap(neighbor, i, j);
                int neighborDistance = LenCount(neighbor, matrix);

                if (neighborDistance < bestNeighborDistance) {
                    bestNeighbor = neighbor;
                    bestNeighborDistance = neighborDistance;
                    improved = true;
                }
            }
        }*/

        // Если нашли улучшение, переходим к лучшему соседу
        if (improved) {
            currentPath = bestNeighbor;
            currentDistance = bestNeighborDistance;
            //PrintVec(bestNeighbor);
            //std::cout << std::endl << "score " << bestNeighborDistance << std::endl;
        }
        else {
            break; // Локальный оптимум достигнут
        }
    }


    return currentPath;
}
std::vector<int> multiStartsteepestHillClimbing(int numRestarts, std::vector<std::vector<double>>& matrix) {
    std::vector<int> bestPath;
    double bestDistance = std::numeric_limits<double>::max();

    // Инициализация генератора случайных чисел
    

    for (int i = 0; i < numRestarts; ++i) {
       
        // Запускаем Hill Climbing
        std::vector<int> localOptimum = hillsteepestClimbingTSP(100, matrix);
        //std::cout <<"\n" << "----------------------------" << std::endl;
        double currentDistance  = LenCount(localOptimum, matrix);

        // Сохраняем лучший результат
        if (currentDistance < bestDistance) {
            bestDistance = currentDistance;
            bestPath = localOptimum;
        }
    }

    return bestPath;
}
std::vector<int> multiStartHillClimbing(int numRestarts, std::vector<std::vector<double>>& matrix) {
    std::vector<int> bestPath;
    double bestDistance = std::numeric_limits<double>::max();

    // Инициализация генератора случайных чисел


    for (int i = 0; i < numRestarts; ++i) {

        // Запускаем Hill Climbing
        std::vector<int> localOptimum = hillClimbingTSP(100, matrix);
        //std::cout <<"\n" << "----------------------------" << std::endl;
        double currentDistance = LenCount(localOptimum, matrix);

        // Сохраняем лучший результат
        if (currentDistance < bestDistance) {
            bestDistance = currentDistance;
            bestPath = localOptimum;
        }
    }

    return bestPath;
}

int main()
{
    size_t n;
    //std::cin >> n;
    std::vector<std::vector<double>> roads = {
        { 0, 14, 16, 41, 37, 39},
        {21,  0, 42, 33, 19, 20},
        { 8, 20,  0, 27, 46, 49},
        {13, 33, 20,  0, 23,  6},
        {12, 32, 43, 24,  0,  4},
        {19, 49, 48, 47, 45,  0}
    };
    
    //FillMatrix(roads); 
    //PrintMatrix(roads);
    //roads={              }
    std::cout << std::endl << "start timer" << std::endl;
    auto start = std::chrono::steady_clock::now();
    std::vector<int> route1 = multiStartHillClimbing(10, roads);
    auto end = std::chrono::steady_clock::now();
    std::cout << std::endl << "end timer" << std::endl;
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << std::endl << duration.count()<< std::endl;
    PrintVec(route1);
    int score = LenCount(route1, roads);
    std::cout << std::endl << score;

    std::cout << std::endl << "start timer" << std::endl;
    start = std::chrono::steady_clock::now();
    route1 = multiStartsteepestHillClimbing(10, roads);
    end = std::chrono::steady_clock::now();
    std::cout << std::endl << "end timer" << std::endl;
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << std::endl << duration.count() << std::endl;
    PrintVec(route1);
    score = LenCount(route1, roads);
    std::cout << std::endl << score;
    //std::cout << "___________________________________" << "\n";
    //std::vector<int> route2 = hillsteepestClimbingTSP(100, roads);
    //PrintVec(route2);
    //score = LenCount(route2, roads);

    //std::cout<< std::endl<< score;

}
 
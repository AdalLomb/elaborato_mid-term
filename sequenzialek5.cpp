#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>
#include <fstream>
#include <sstream>

using namespace std;

// Calcola la distanza Euclidea tra due punti
double distance(vector<double>& a, vector<double>& b) {
    double dist = 0.0;
    for (int i = 0; i < a.size(); i++) {
        dist += pow(a[i] - b[i], 2);
    }
    return sqrt(dist);
}

// Esegue l'algoritmo K-means clustering
void kMeansClustering(vector<vector<double>>& data, int k, vector<vector<double>>& centroids, vector<int>& clusters) {
    int n = data.size();
    int dim = data[0].size();
    
    // Inizializza i centroidi in modo casuale
    mt19937 gen(time(0));
    uniform_int_distribution<> dist(0, n-1);
    for (int i = 0; i < k; i++) {
        int idx = dist(gen);
        centroids[i] = data[idx];
    }
    
    // Assegna ogni punto al cluster più vicino
    for (int i = 0; i < n; i++) {
        double minDist = 5;
        int minIdx = 0;
        for (int j = 0; j < k; j++) {
            double dist = distance(data[i], centroids[j]);
            if (dist < minDist) {
                minDist = dist;
                minIdx = j;
            }
        }
        clusters[i] = minIdx;
    }
    
    // Ripete fino a convergenza
    bool changed = true;
    int iter = 0;
    while (changed && iter < 100) {
        // Calcola i nuovi centroidi
        vector<int> counts(k, 0);
        vector<vector<double>> newCentroids(k, vector<double>(dim, 0.0));
        for (int i = 0; i < n; i++) {
            int cluster = clusters[i];
            counts[cluster]++;
            for (int j = 0; j < dim; j++) {
                newCentroids[cluster][j] += data[i][j];
            }
        }
        for (int i = 0; i < k; i++) {
            if (counts[i] > 0) {
                for (int j = 0; j < dim; j++) {
                    newCentroids[i][j] /= counts[i];
                }
            }
        }
        
        // Assegna ogni punto al cluster più vicino
        changed = false;
        for (int i = 0; i < n; i++) {
            double minDist = 5;
            int minIdx = 0;
            for (int j = 0; j < k; j++) {
                double dist = distance(data[i], newCentroids[j]);
                if (dist < minDist) {
                    minDist = dist;
                    minIdx = j;
                }
            }
            if (minIdx != clusters[i]) {
                changed = true;
                clusters[i] = minIdx;
            }
        }
        centroids = newCentroids;
        iter++;
    }
}


	int main() {
    // Apre il file di input
    std::ifstream inputFile("sample_1000000.txt");

    // Crea il vector che conterrà i dati
    std::vector<std::vector<double>> data;

    // Legge il contenuto del file riga per riga
    std::string line;
    std::getline(inputFile, line);

    // Crea uno stringstream per dividere la riga in sottostringhe
    std::stringstream ss(line);

    // Legge i valori tra parentesi graffe e li aggiunge al vector "data"
    char c;
    while (ss >> c) {
        if (c == '{') {
            std::vector<double> point;
            double x, y;
            ss >> x >> c >> y >> c;
            point.push_back(x);
            point.push_back(y);
            data.push_back(point);
        }
    }

    // Chiude il file di input
    inputFile.close();
	
		int k = 5; 
		// Centroidi
		vector<vector<double>> centroids(k, vector<double>(2));
		// Cluster	
		vector<int> clusters(data.size()); 
		
		// Misura il tempo di elaborazione
		clock_t start = clock();
		
		// Esegue l'algoritmo K-means clustering
		kMeansClustering(data, k, centroids, clusters);
		
		// Calcola il tempo di elaborazione
		clock_t end = clock();
		double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
		
		 cout << "Tempo di elaborazione sequenziale: " << elapsed_secs << " secondi" << endl;
		 
		return 0;
	}

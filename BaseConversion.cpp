#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>

int main() {
    double N;
    double B_input;
    int P;

    std::cout << "Enter N: "; std::cin >> N;
    std::cout << "Enter B: "; std::cin >> B_input;
    while (B_input <= 1 || B_input > 10 || B_input != std::floor(B_input)) {
        std::cout << "B cannot be 1, less than 1, or above 10, enter new B: "; std::cin >> B_input;
    }
    std::cout << "Enter P: "; std::cin >> P;

    int B = static_cast<int>(B_input);
    double I = std::floor(N);
    std::string convertedInt= "";
    long long tempInt = static_cast<long long>(I);

    if (tempInt == 0) {
        convertedInt = "0";
    } else {
        while (tempInt > 0) {
            convertedInt += "0123456789"[tempInt%B];
            tempInt /= B;
        }
    std::reverse(convertedInt.begin(), convertedInt.end());
    }

    std::string convertedFrac = "";
    double F = N - I;
    if (F == 0){
        convertedFrac = "0";
    } else {
        while (F > 1e-9) {
            double curr = B * F;
            int carryWI = std::floor(curr);
            convertedFrac += "0123456789"[carryWI];
            F = curr - carryWI;
        }
    }

    std::string total = convertedInt + '.' + convertedFrac;

    double returnedInt = 0;
    double returnedFrac = 0;
    double Bpower = B;

    std::string truncatedFrac = "";
    std::string truncatedInt = "";
    int digitCount = 0;
    bool pass = false;
    for (char c : total) {
        if (c == '.') {
            pass = true;
            continue;
        }

        digitCount++;

        if (!pass) {
            truncatedInt += c;
        } else {
            truncatedFrac += c;
        }

        if (digitCount == P) {
            break;
        }
    }

    for (int i=0; i < truncatedInt.length(); i++) {
        int digit = truncatedInt[i] - '0';
        returnedInt = (returnedInt*B) + digit;
    }

    for (int i=0; i < truncatedFrac.length(); i++) {
        double digit = truncatedFrac[i] - '0';
        returnedFrac += digit / Bpower;
        Bpower = Bpower * B;
    }

    double returned = returnedInt + returnedFrac;
    double epsilon = (std::abs(N - returned)) / (N == 0 ? 1 : N);

    std::cout << "a) Floating Point Representation: " << total << std::endl;
    std::cout << "b)i) Returned Number: " << returned << std::endl;
    std::cout << "b)ii) Relative Error: " << epsilon << std::endl;

    return 0;
}
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

using namespace std;
int main() {
    string buf;
    cin >> buf;
    vector<string> vs;
    for (int i = 0; i < buf.length(); ++i) {
        buf = buf[buf.length()-1] + buf.substr(0, buf.length() - 1);
        vs.push_back(buf);
    }
    sort(vs.begin(), vs.end());
    for (int i = 0; i < vs.size(); ++i) {
        cout << vs[i] << endl;
        cerr << vs[i][vs[i].length() - 1];
    }
    cerr << endl;
    return 0;
}

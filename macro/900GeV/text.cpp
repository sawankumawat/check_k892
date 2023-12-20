#include<iostream>
using namespace std;
int main(){
    int a;

    const char* hi[2] = {"hello1", "hello2"};
    if (hi[0] == "hello1")
    {
         a = 69;
    }
    cout << a<<endl;
    

    return 0;
}
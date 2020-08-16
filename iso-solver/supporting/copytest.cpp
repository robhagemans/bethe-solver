#include <iostream>
using std::cout;
using std::endl;

class CopyTest {
public:
    CopyTest();
    ~CopyTest();
    CopyTest(const CopyTest& orig);
    CopyTest& operator=(const CopyTest& orig);
    int z;
};


CopyTest::CopyTest()
{
    cout<<" init"<<endl;
    z=1;
}


CopyTest::~CopyTest()
{
    cout<<" destroy"<<endl;
}

CopyTest::CopyTest(const CopyTest& orig) 
{
    cout<<" copy "<<endl;
    z = orig.z;
}


CopyTest& CopyTest::operator=(const CopyTest& orig) 
{
    cout<<" assign "<<endl;
    if (&orig==this) return (*this);
    z= orig.z;
    return (*this);
}

CopyTest TestReturn()
{
    CopyTest a;
    return a;
}


CopyTest TestReturn2()
{
    return CopyTest();
}

int main()
{
    cout<<"on return: "<<endl;
    CopyTest b = TestReturn();

    cout<<"on return 2: "<<endl;
    CopyTest d = TestReturn2();

    cout<<"explicit copy: "<<endl;
    CopyTest c = b;
    
    
    cout<<"explicit assign: "<<endl;
    b = c;

    cout<<"end: "<<endl;
    return 0;
}




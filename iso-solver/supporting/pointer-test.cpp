
void a(int* test)
{
    delete test;
}


void b(int*& test)
{
    delete test;
    test=0;
}

void c(int** test)
{
    delete (*test);
    *test=0;
}



int main()
{
    a(new int(3));
    c(&(new int(3)));
    
    int x=5;
    c(&&x);
    
    return 0;
}


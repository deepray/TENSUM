#include "safe_logic.h"

// SafeLess
bool SafeLess( int x, int y )
{
    return (x < y);
}

bool SafeLess( int x, unsigned int y )
{
    return (x < 0) || (static_cast<unsigned int>(x) < y);
}

bool SafeLess( int x, unsigned long int y )
{
    return (x < 0) || (static_cast<unsigned long int>(x) < y);
}

bool SafeLess( unsigned int x, int y )
{
    return (y >= 0) && (x < static_cast<unsigned int>(y));
}

bool SafeLess( unsigned int x, unsigned int y )
{
    return (x < y);
}

bool SafeLess( unsigned int x, unsigned long int y )
{
    return (static_cast<unsigned long int>(x) < y);
}

bool SafeLess( unsigned long int x, int y )
{
    return (y >= 0) && (x < static_cast<unsigned long int>(y));
}

bool SafeLess( unsigned long int x, unsigned int y )
{
    return ( x < static_cast<unsigned long int>(y));
}

bool SafeLess( unsigned long int x, unsigned long int y )
{
    return (x < y);
}

// SafeLessEq
bool SafeLessEq( int x, int y )
{
    return (x <= y);
}

bool SafeLessEq( int x, unsigned int y )
{
    return (x < 0) || (static_cast<unsigned int>(x) <= y);
}

bool SafeLessEq( int x, unsigned long int y )
{
    return (x < 0) || (static_cast<unsigned long int>(x) <= y);
}

bool SafeLessEq( unsigned int x, int y )
{
    return (y >= 0) && (x <= static_cast<unsigned int>(y));
}

bool SafeLessEq( unsigned int x, unsigned int y )
{
    return (x <= y);
}

bool SafeLessEq( unsigned int x, unsigned long int y )
{
    return (static_cast<unsigned long int>(x) <= y);
}

bool SafeLessEq( unsigned long int x, int y )
{
    return (y >= 0) && (x <= static_cast<unsigned long int>(y));
}

bool SafeLessEq( unsigned long int x, unsigned int y )
{
    return ( x <= static_cast<unsigned long int>(y));
}

bool SafeLessEq( unsigned long int x, unsigned long int y )
{
    return (x <= y);
}



// SafeMore
bool SafeMore( int x, int y )
{
    return SafeLess(y,x);
}

bool SafeMore( int x, unsigned int y )
{
    return SafeLess(y,x);
}

bool SafeMore( int x, unsigned long int y )
{
    return SafeLess(y,x);
}

bool SafeMore( unsigned int x, int y )
{
    return SafeLess(y,x);
}

bool SafeMore( unsigned int x, unsigned int y )
{
    return SafeLess(y,x);
}

bool SafeMore( unsigned int x, unsigned long int y )
{
    return SafeLess(y,x);
}

bool SafeMore( unsigned long int x, int y )
{
    return SafeLess(y,x);
}

bool SafeMore( unsigned long int x, unsigned int y )
{
    return SafeLess(y,x);
}

bool SafeMore( unsigned long int x, unsigned long int y )
{
    return SafeLess(y,x);
}



// SafeMoreEq
bool SafeMoreEq( int x, int y )
{
    return SafeLessEq(y,x);
}

bool SafeMoreEq( int x, unsigned int y )
{
    return SafeLessEq(y,x);
}

bool SafeMoreEq( int x, unsigned long int y )
{
    return SafeLessEq(y,x);
}

bool SafeMoreEq( unsigned int x, int y )
{
    return SafeLessEq(y,x);
}

bool SafeMoreEq( unsigned int x, unsigned int y )
{
    return SafeLessEq(y,x);
}

bool SafeMoreEq( unsigned int x, unsigned long int y )
{
    return SafeLessEq(y,x);
}

bool SafeMoreEq( unsigned long int x, int y )
{
    return SafeLessEq(y,x);
}

bool SafeMoreEq( unsigned long int x, unsigned int y )
{
    return SafeLessEq(y,x);
}

bool SafeMoreEq( unsigned long int x, unsigned long int y )
{
    return SafeLessEq(y,x);
}


// SafeEq
bool SafeEq( int x, int y )
{
    return (x == y);
}

bool SafeEq( int x, unsigned int y )
{
    return (x >= 0) && (static_cast<unsigned int>(x) == y);
}

bool SafeEq( int x, unsigned long int y )
{
    return (x >= 0) || (static_cast<unsigned long int>(x) == y);
}

bool SafeEq( unsigned int x, int y )
{
    return (y >= 0) && (x == static_cast<unsigned int>(y));
}

bool SafeEq( unsigned int x, unsigned int y )
{
    return (x == y);
}

bool SafeEq( unsigned int x, unsigned long int y )
{
    return (static_cast<unsigned long int>(x) == y);
}

bool SafeEq( unsigned long int x, int y )
{
    return (y >= 0) && (x == static_cast<unsigned long int>(y));
}

bool SafeEq( unsigned long int x, unsigned int y )
{
    return ( x == static_cast<unsigned long int>(y));
}

bool SafeEq( unsigned long int x, unsigned long int y )
{
    return (x == y);
}


// SafeNeq
bool SafeNeq( int x, int y )
{
    return (!SafeEq(x,y));
}

bool SafeNeq( int x, unsigned int y )
{
    return (!SafeEq(x,y));
}

bool SafeNeq( int x, unsigned long int y )
{
    return (!SafeEq(x,y));
}

bool SafeNeq( unsigned int x, int y )
{
    return (!SafeEq(x,y));
}

bool SafeNeq( unsigned int x, unsigned int y )
{
    return (!SafeEq(x,y));
}

bool SafeNeq( unsigned int x, unsigned long int y )
{
    return (!SafeEq(x,y));
}

bool SafeNeq( unsigned long int x, int y )
{
    return (!SafeEq(x,y));
}

bool SafeNeq( unsigned long int x, unsigned int y )
{
    return (!SafeEq(x,y));
}

bool SafeNeq( unsigned long int x, unsigned long int y )
{
    return (!SafeEq(x,y));
}
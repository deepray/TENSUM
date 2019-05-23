#ifndef __SAFE_LOGIC_H__
#define __SAFE_LOGIC_H__

bool SafeLess( int x, int y );
bool SafeLess( int x, unsigned int y );
bool SafeLess( int x, unsigned long int y );
bool SafeLess( unsigned int x, int y );
bool SafeLess( unsigned int x, unsigned int y );
bool SafeLess( unsigned int x, unsigned long int y );
bool SafeLess( unsigned long int x, int y );
bool SafeLess( unsigned long int x, unsigned int y );
bool SafeLess( unsigned long int x, unsigned long int y );

bool SafeLessEq( int x, int y );
bool SafeLessEq( int x, unsigned int y );
bool SafeLessEq( int x, unsigned long int y );
bool SafeLessEq( unsigned int x, int y );
bool SafeLessEq( unsigned int x, unsigned int y );
bool SafeLessEq( unsigned int x, unsigned long int y );
bool SafeLessEq( unsigned long int x, int y );
bool SafeLessEq( unsigned long int x, unsigned int y );
bool SafeLessEq( unsigned long int x, unsigned long int y );

bool SafeMore( int x, int y );
bool SafeMore( int x, unsigned int y );
bool SafeMore( int x, unsigned long int y );
bool SafeMore( unsigned int x, int y );
bool SafeMore( unsigned int x, unsigned int y );
bool SafeMore( unsigned int x, unsigned long int y );
bool SafeMore( unsigned long int x, int y );
bool SafeMore( unsigned long int x, unsigned int y );
bool SafeMore( unsigned long int x, unsigned long int y );

bool SafeMoreEq( int x, int y );
bool SafeMoreEq( int x, unsigned int y );
bool SafeMoreEq( int x, unsigned long int y );
bool SafeMoreEq( unsigned int x, int y );
bool SafeMoreEq( unsigned int x, unsigned int y );
bool SafeMoreEq( unsigned int x, unsigned long int y );
bool SafeMoreEq( unsigned long int x, int y );
bool SafeMoreEq( unsigned long int x, unsigned int y );
bool SafeMoreEq( unsigned long int x, unsigned long int y );


bool SafeEq( int x, int y );
bool SafeEq( int x, unsigned int y );
bool SafeEq( int x, unsigned long int y );
bool SafeEq( unsigned int x, int y );
bool SafeEq( unsigned int x, unsigned int y );
bool SafeEq( unsigned int x, unsigned long int y );
bool SafeEq( unsigned long int x, int y );
bool SafeEq( unsigned long int x, unsigned int y );
bool SafeEq( unsigned long int x, unsigned long int y );

bool SafeNeq( int x, int y );
bool SafeNeq( int x, unsigned int y );
bool SafeNeq( int x, unsigned long int y );
bool SafeNeq( unsigned int x, int y );
bool SafeNeq( unsigned int x, unsigned int y );
bool SafeNeq( unsigned int x, unsigned long int y );
bool SafeNeq( unsigned long int x, int y );
bool SafeNeq( unsigned long int x, unsigned int y );
bool SafeNeq( unsigned long int x, unsigned long int y );


#endif
#ifndef mymacros_h
#define mymacros_h

//Shorthands

#define For(i,N) for(int i=0;i<N;i++)

#define ARRAY_COUNT( array ) (sizeof( array ) / (sizeof( array[0] ) * (sizeof( array ) != sizeof(void*) || sizeof( array[0] ) <= sizeof(void*))))

//#undef For(i,N) for(int i=0;i<N;i++)

//#define For(i,a,N) for(int i=a;i<N;i++)
//#undef For(i,a,N) for(int i=a;i<N;i++)

//#define For(i,a,b,N) for(int i=a;i<N;i+=b)
//#undef For(i,a,b,N) for(int i=a;i<N;i+=b)


#endif

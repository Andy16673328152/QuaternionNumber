# FAQ

Q:Why not GLM/Eigen?

A:GLM/Eigen don't provide transcendental functions.If you want to use them,choose me.

Q:Code by whom?

A:Me.AI wouldn't code N releases,fix bug of cosh,and CI test.

Q:Why operator>> use string parsing,not regex?

A:Regex cannot parse 1+-i.

Q:Is it possible that _mm_load_ps collapses?

A:No.The struct has alignas(16),in most systems,std::max_alignof=16

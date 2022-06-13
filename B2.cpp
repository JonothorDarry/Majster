#include <stdio.h>
#include <vector>
#include <cmath>
#include "common.hpp"
using ll=long long;
using namespace std;
const int D=50001, Z=701; //Metaconst - 701 - remove!!!

ll dp[Z][D];
void generate_dp(ll limit, vector <ll> &prime){
	for (ll j=0; j<D; j++) dp[0][j] = j;
	for (ll i=1; i<=limit; i++){
		for (ll j=0; j<D; j++) dp[i][j] = dp[i-1][j] - dp[i-1][j/prime[i]];
	}
}

ll calc_phi(ll x, ll limit, vector <ll> &prime){
	if (x<D) return dp[limit][x];
	if (limit == 0) return x;
	return calc_phi(x, limit-1, prime) - calc_phi(x/prime[limit], limit-1, prime);
}

ll phi_part(ll x, ll limit, vector<ll> &prime){
	return calc_phi(x, limit, prime);
}

vector <ll> prime;
ll prime_counter(ll x){
	if (x==1) return 0;

	ll limiting_prime = get_prime_above(prime, ll(ceil(pow(x, 1/3.)))); //to binary searcher
	return limiting_prime - 1 + phi_part(x, limiting_prime, prime) - part_two(x, limiting_prime+1, prime) - part_three(x, limiting_prime+1, prime);
}

void shared_part(ll maximal_x){
	ll last_prime = ll(ceil(pow(maximal_x, 2/3.)))+1;
	get_primes(last_prime, prime);
	generate_dp(Z-1, prime);
}

int main(){
	shared_part(MAXIM);
	int t; scanf ("%d", &t);
	while(t--){
		ll x; scanf ("%lld", &x);
		printf ("%lld\n", prime_counter(x));
	}
return 0;}

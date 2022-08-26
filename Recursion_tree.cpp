#include <stdio.h>
#include <vector>
#include <cmath>
#include "common.hpp"
using ll=long long;
using namespace std;
const int D=100001;

vector <vector <pair<ll,ll>> > layer(D);
ll phi_part(ll x, ll limit, vector<ll> &prime){
	for (int i=limit; i>=0; i--) layer[i].clear();

	layer[limit].push_back({x, 0});
	for (int i=limit; i>0; i--){
		int ln = layer[i].size();
		int ite_1 = 0, ite_2 = 0;

		ll last_elem = -1;
		while (ite_1 < ln || ite_2 < ln){
			ll elem_1=Inf, elem_2=Inf;
			if (ite_1 < ln) elem_1 = layer[i][ite_1].first;
			if (ite_2 < ln) elem_2 = layer[i][ite_2].first/prime[i];

			ll next_elem = min(elem_1, elem_2);
			if (next_elem != last_elem){
				layer[i-1].push_back({next_elem, 0});
				last_elem = next_elem;
			}

			if (elem_1 == elem_2) ite_1++, ite_2++;
			else if (elem_1 < elem_2) ite_1++;
			else ite_2++;
		}
	}

	for (int j=0; j<layer[0].size(); j++) layer[0][j].second = layer[0][j].first;
	for (int i=1; i<=limit; i++){
		for (int j=0, k=0; k<layer[i].size(); k++){
			while (layer[i][k].first != layer[i-1][j].first) j++;
			layer[i][k].second += layer[i-1][j].second;
		}

		for (int j=0, k=0; k<layer[i].size(); k++){
			while (layer[i][k].first/prime[i] != layer[i-1][j].first) j++;
			layer[i][k].second -= layer[i-1][j].second;
		}
	}

	return layer[limit][0].second;
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
}

int main(){
	shared_part(MAXIM);
	int t; scanf ("%d", &t);
	while(t--){
		ll x; scanf ("%lld", &x);
		printf ("%lld\n", prime_counter(x));
	}
return 0;}

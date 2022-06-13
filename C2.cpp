#include <stdio.h>
#include <vector>
#include <cmath>
#include "common.hpp"
using ll=long long;
using namespace std;
const int D=100001, Cv=1<<25; //Independent const!

int val[Cv<<1];
void initialize_segment_tree(ll limit){
	val[Cv] = 0;
	for (int i=1; i<=limit; i+=2) val[Cv + i] = 1;
	for (int i=Cv-1; i>0; i--) val[i] = val[i<<1] + val[(i<<1)+1];
}

void update_segment_tree(ll basis, ll limit){
	for (ll value=basis; value <= limit; value+=basis){
		if (val[Cv + value] == 0) continue;
		for (ll x=Cv+value; x>0; x>>=1) val[x]--;
	}
}

ll query_segment_tree(ll l, ll r){
	ll res = val[l+Cv];
	if (l != r) res += val[r+Cv];

	for (l+=Cv,r+=Cv; l>0; l>>=1,r>>=1){
		if (l%2==0 && r-l>1) res += val[l+1];
		if (r%2==1 && r-l>1) res += val[r-1];
	}
	return res;
}

vector <vector <pair<ll,ll>> > layer(D);
ll phi_part(ll x, ll limit, vector<ll> &prime){
	for (int i=limit; i>=0; i--) layer[i].clear();
	ll proper_limit = x/prime[limit];

	layer[limit].push_back({x, 0});
	for (int i=limit; i>0; i--){
		int ln = layer[i].size();
		int ite_1 = 0, ite_2 = 0;

		ll last_elem = -1;
		for (; layer[i][ite_2].first/prime[i] <= proper_limit && ite_2 < ln; ite_2++) ;

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

	initialize_segment_tree(proper_limit);
	for (int j=0; j<layer[0].size(); j++) layer[0][j].second = layer[0][j].first;
	for (int j=0; j<layer[1].size(); j++) layer[1][j].second = layer[1][j].first/2 + layer[1][j].first%2;

	for (int i=2; i<=limit; i++){
		for (int j=0, k=0; k<layer[i].size(); k++){
			while (layer[i][k].first != layer[i-1][j].first) j++;
			layer[i][k].second += layer[i-1][j].second;
		}

		for (int j=0, k=0; k<layer[i].size(); k++){
			if (layer[i][k].first/prime[i] <= proper_limit) layer[i][k].second -= query_segment_tree(0, layer[i][k].first/prime[i]);
			else{
				while (layer[i][k].first/prime[i] != layer[i-1][j].first) j++;
				layer[i][k].second -= layer[i-1][j].second;
			}
		}
		update_segment_tree(prime[i], proper_limit);
	}

	return layer[limit][0].second;
}

vector <ll> prime;
ll prime_counter(ll x){
	if (x==1) return 0;

	ll limiting_prime = get_prime_above(prime, ll(ceil(pow(x, 1/3.)))); //to binary searcher
	return limiting_prime - 1 + phi_part(x, limiting_prime, prime) - part_two(x, limiting_prime+1, prime);
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

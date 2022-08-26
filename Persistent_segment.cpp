#include <stdio.h>
#include <vector>
#include <cmath>
#include "common.hpp"
using ll=long long;
using namespace std;
const int BCv=25, Cv=1<<BCv;

struct segtree_vertex{
	int value;
	int left;
	int right;
	int index;
};

vector <vector <segtree_vertex> > val(Cv<<1);
void initialize_segment_tree(ll limit){
	for (int i=0; i<Cv; i++) val[Cv + i].reserve(2);

	val[Cv].push_back({0, -1, -1, 0});
	for (int i=1; i<=limit; i++) val[Cv + i].push_back({1, -1, -1, 0});
	for (int i=limit+1; i<Cv; i++) val[Cv + i].push_back({0, -1, -1, 0});
}

void update_segment_tree(int basis, ll limit){
	for (ll value=basis; value <= limit; value+=basis){
		if (val[Cv + value].back().value == 0) continue;
		val[value+Cv].push_back({0, -1, -1, basis});
	}
}

void after_all_updates(){
	for (ll x=Cv-1; x>0; x--){
		int nexte = (x<<1);
		int l_length = val[nexte].size(), r_length = val[nexte+1].size();

		val[x].reserve(max(l_length, r_length)+1);
		for (int l=0,r=0; l<l_length || r<r_length; ){
			int index=-1;

			if (l < l_length && r < r_length && val[nexte+1][r].index == val[nexte][l].index) index = val[nexte][l].index, l++, r++;
			else if (l == l_length || (r < r_length && val[nexte+1][r].index < val[nexte][l].index)) index = val[nexte+1][r].index, r++;
			else index = val[nexte][l].index, l++;

			val[x].push_back({val[nexte][l-1].value + val[nexte+1][r-1].value, l-1, r-1, index});
		}
	}
}

ll query_segment_tree(ll l, ll r, int prime_nr){
	pair<int,int> node_left={1, prime_nr}, node_right = {1, prime_nr}; //(x,y)
	int res=0, bit = Cv>>1, l_root, r_root;

	int l_x, r_x;
	while (true){
		l_x = node_left.first, r_x = node_right.first;
		if (l_x%2 == 0 && r_x-l_x > 1) res += val[l_x+1][val[l_x>>1][l_root].right].value;
		if (r_x%2 == 1 && r_x-l_x > 1) res += val[r_x-1][val[r_x>>1][r_root].left].value;

		if (bit <= 0) break;
		l_root = node_left.second;
		r_root = node_right.second;

		if ((l&bit) == 0) node_left = {node_left.first<<1, val[node_left.first][node_left.second].left};
		else node_left = {(node_left.first<<1)+1, val[node_left.first][node_left.second].right};

		if ((r&bit) == 0) node_right = {node_right.first<<1, val[node_right.first][node_right.second].left};
		else node_right = {(node_right.first<<1)+1, val[node_right.first][node_right.second].right};

		bit>>=1;
	}
	res += val[l_x][node_left.second].value;
	if (l_x != r_x) res += val[r_x][node_right.second].value;

	return res;
}

void make_segment_tree(ll proper_limit, ll limit, vector <ll> &prime){
	initialize_segment_tree(proper_limit);
	for (int i=1; i<=limit; i++) update_segment_tree(prime[i], proper_limit);
	after_all_updates();
}

ll rec_phi(ll x, ll a, ll proper_limit, vector<ll> &prime){
	if (a == 0) return x;
	if (x <= proper_limit) return query_segment_tree(0, x, a);
	return rec_phi(x, a-1, proper_limit, prime) - rec_phi(x/prime[a], a-1, proper_limit, prime);

}

ll phi_part(ll x, ll limit, vector<ll> &prime){
	ll proper_limit = x/prime[limit];
	return rec_phi(x, limit, proper_limit, prime);
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
	ll limiting_prime = get_prime_above(prime, ll(ceil(pow(maximal_x, 1/3.)))); //to binary searcher

	ll proper_limit = maximal_x/prime[limiting_prime];
	make_segment_tree(proper_limit, limiting_prime, prime);
}

int main(){
	shared_part(MAXIM);
	int t; scanf ("%d", &t);
	while(t--){
		ll x; scanf ("%lld", &x);
		printf ("%lld\n", prime_counter(x));
	}
return 0;}

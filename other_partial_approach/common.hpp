#include <vector>
#include <algorithm>
using ll=long long;
using namespace std;
const ll Inf=1000000000000000001ll, MAXIM=100000000001ll;
const int PRIME_CACHE=50001;

bool is_prime[PRIME_CACHE];
void get_primes(ll range, vector<ll> &prime){
	ll basis = ll(ceil(pow(range, 1/2.)));

	for (ll i=2; i<basis; i++) is_prime[i] = true;
	for (ll i=2; i*i<basis; i++){
		if (!is_prime[i]) continue;
		for (ll j=i*i; j<basis; j+=i) is_prime[j] = false;
	}

	prime.reserve(ceil(range/log2(range)));
	prime.push_back(-1);
	for (ll i=0; i < basis; i++){
	       if (is_prime[i]) prime.push_back(i);
	}
	ll last_primal_prime = prime.size();

	for (ll k=1; k*basis < range; k++){
		ll state = k*basis;
		for (ll i=0; i < basis; i++) is_prime[i] = true;

		for (ll i=1; i<last_primal_prime; i++){
			ll mapped_basis = prime[i] - ((state%prime[i] == 0)?prime[i]:(state % prime[i]));
			for (ll j=mapped_basis; j<basis; j+=prime[i]) is_prime[j] = false;
		}

		for (ll i=0; i < basis; i++){
		       if (is_prime[i] && state+i < range) prime.push_back(state+i);
		}
	}
}

ll part_two(ll x, ll limit, vector<ll> &prime){
	ll j=prime.size()-1;
	ll res=0;
	for (ll i=limit; i <= j; i++){
		while (prime[j]*prime[i] > x) j--;
		if (j<i) break;
		res += j - i + 1;
	}
	return res;
}

ll part_three(ll x, ll limit, vector<ll> &prime){
	ll res=0, last_k=prime.size()-1;
	for (ll i=limit; prime[i] * prime[i] * prime[i] <= x; i++){
		ll k = last_k;
		while (prime[i]*prime[i]*prime[k] > x) k--;
		last_k = k;

		for (ll j=i; j <= k; j++){
			while (prime[j]*prime[i]*prime[k] > x) k--;
			if (k<j) break;
			res += k - j + 1;
		}
	}
	return res;
}

//lowest ge
ll get_prime_above(vector <ll> &primes, ll limit){
	return std::lower_bound(primes.begin(), primes.end(), limit) - primes.begin();
}

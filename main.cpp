#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <set>
#include <utility>

struct Sekwencja {
    std::string id;
    std::string sekwencja;
};

struct Jakosc {
    std::string id;
    std::vector<int> oceny_jakosci;
};

struct WierzcholekGrafu {
    std::string idSekwencji;
    std::string podciag;
    int pozycja; // Pozycja podciągu w sekwencji
    bool operator<(const WierzcholekGrafu& other) const {
        if (idSekwencji != other.idSekwencji) return idSekwencji < other.idSekwencji;
        if (podciag != other.podciag) return podciag < other.podciag;
        return pozycja < other.pozycja;
    }
};

struct PrzetworzonaSekwencja {
    std::string id;
    std::string sekwencja;
    std::vector<int> pozycje;
    std::vector<int> usunietePozycje;
    std::vector<int> pozycjeOryginalne;
    std::map<int, int> usunieteNukleotydy;
};

struct Krawedz {
    WierzcholekGrafu wierzcholek1;
    WierzcholekGrafu wierzcholek2;

    Krawedz(const WierzcholekGrafu& w1, const WierzcholekGrafu& w2) : wierzcholek1(w1), wierzcholek2(w2) {}

    bool operator<(const Krawedz& other) const {
        if (wierzcholek1 < other.wierzcholek1) return true;
        if (other.wierzcholek1 < wierzcholek1) return false;
        return wierzcholek2 < other.wierzcholek2;
    }
};

std::map<std::string, Sekwencja> wczytajFasta() {
    std::map<std::string, Sekwencja> sekwencje;
    std::ifstream file("fasta.txt");
    std::string linia;
    Sekwencja aktualnaSekwencja;

    while (getline(file, linia)) {
        if (linia[0] == '>') {
            if (!aktualnaSekwencja.id.empty()) {
                sekwencje[aktualnaSekwencja.id] = aktualnaSekwencja;
                aktualnaSekwencja.sekwencja.clear();
            }
            aktualnaSekwencja.id = linia.substr(1, linia.find(" ") - 1);
        } else {
            aktualnaSekwencja.sekwencja += linia;
        }
    }
    if (!aktualnaSekwencja.id.empty()) {
        sekwencje[aktualnaSekwencja.id] = aktualnaSekwencja;
    }

    return sekwencje;
}

std::map<std::string, Jakosc> wczytajJakosc() {
    std::map<std::string, Jakosc> jakosci;
    std::ifstream file("wiarygodnosc.txt");
    std::string linia;
    Jakosc aktualnaJakosc;

    while (getline(file, linia)) {
        if (linia[0] == '>') {
            if (!aktualnaJakosc.id.empty()) {
                jakosci[aktualnaJakosc.id] = aktualnaJakosc;
                aktualnaJakosc.oceny_jakosci.clear();
            }
            aktualnaJakosc.id = linia.substr(1, linia.find(" ") - 1);
        } else {
            std::istringstream iss(linia);
            int ocena;
            while (iss >> ocena) {
                aktualnaJakosc.oceny_jakosci.push_back(ocena);
            }
        }
    }
    if (!aktualnaJakosc.id.empty()) {
        jakosci[aktualnaJakosc.id] = aktualnaJakosc;
    }

    return jakosci;
}

std::map<std::string, PrzetworzonaSekwencja> usunNiskiejJakosciNukleotydy(
        const std::map<std::string, Sekwencja>& sekwencje,
        const std::map<std::string, Jakosc>& jakosci,
        int progJakosci) {

    std::map<std::string, PrzetworzonaSekwencja> przetworzoneSekwencje;

    for (const auto& para : sekwencje) {
        const auto& id = para.first;
        const auto& sekwencja = para.second.sekwencja;
        const auto& oceny_jakosci = jakosci.at(id).oceny_jakosci;

        PrzetworzonaSekwencja przetworzona;
        przetworzona.id = id;

        for (size_t i = 0; i < sekwencja.size(); ++i) {
            if (oceny_jakosci[i] >= progJakosci) {
                przetworzona.sekwencja.push_back(sekwencja[i]);
                przetworzona.pozycje.push_back(i + 1);
                przetworzona.pozycjeOryginalne.push_back(i + 1);
            } else {
                przetworzona.usunieteNukleotydy[i + 1] = oceny_jakosci[i];
            }
        }

        przetworzoneSekwencje[id] = przetworzona;
    }

    return przetworzoneSekwencje;
}


using Graf = std::vector<WierzcholekGrafu>;

Graf utworzGraf(const std::map<std::string, PrzetworzonaSekwencja>& sekwencje, int dlugoscPodciagu) {
    Graf graf;

    for (const auto& para : sekwencje) {
        const auto& idSekwencji = para.first;
        const auto& sekwencja = para.second.sekwencja;
        const auto& pozycjeOryginalne = para.second.pozycjeOryginalne;

        for (size_t i = 0; i <= sekwencja.size() - dlugoscPodciagu; ++i) {
            WierzcholekGrafu wierzcholek;
            wierzcholek.idSekwencji = idSekwencji;
            wierzcholek.podciag = sekwencja.substr(i, dlugoscPodciagu);
            wierzcholek.pozycja = pozycjeOryginalne[i];
            graf.push_back(wierzcholek);
        }
    }

    return graf;
}

using ZbiorKrawedzi = std::set<Krawedz>;

ZbiorKrawedzi polaczWierzcholki(const Graf& graf, int dlugoscPodciagu) {
    ZbiorKrawedzi krawedzie;

    for (size_t i = 0; i < graf.size(); ++i) {
        for (size_t j = i + 1; j < graf.size(); ++j) {
            if (graf[i].idSekwencji != graf[j].idSekwencji &&
                graf[i].podciag == graf[j].podciag &&
                std::abs(graf[i].pozycja - graf[j].pozycja) <= 10 * dlugoscPodciagu) {
                Krawedz krawedz(graf[i], graf[j]);
                krawedzie.insert(krawedz);
            }
        }
    }

    return krawedzie;
}


void wyswietlWierzcholkiGrafu(const Graf& graf) {
    for (const auto& wierzcholek : graf) {
        std::cout << "Wierzchołek-Podciąg: " << wierzcholek.podciag << "\n";
        std::cout << "  Sekwencja ID: " << wierzcholek.idSekwencji << ", Pozycja: " << wierzcholek.pozycja << "\n";
        std::cout << "\n";
    }
}
void znajdzStruktureGwiazdy(const Graf& graf, const ZbiorKrawedzi& krawedzie) {
    for (const auto& wierzcholek : graf) {
        std::set<std::string> unikalneSekwencje;
        std::vector<WierzcholekGrafu> sasiedzi;

        for (const auto& krawedz : krawedzie) {
            if (krawedz.wierzcholek1.idSekwencji == wierzcholek.idSekwencji && krawedz.wierzcholek1.podciag == wierzcholek.podciag) {
                sasiedzi.push_back(krawedz.wierzcholek2);
                unikalneSekwencje.insert(krawedz.wierzcholek2.idSekwencji);
            } else if (krawedz.wierzcholek2.idSekwencji == wierzcholek.idSekwencji && krawedz.wierzcholek2.podciag == wierzcholek.podciag) {
                sasiedzi.push_back(krawedz.wierzcholek1);
                unikalneSekwencje.insert(krawedz.wierzcholek1.idSekwencji);
            }
        }

        if (sasiedzi.size() >= 4 && unikalneSekwencje.size() >= 4) {
            std::cout << "Znaleziono strukturę typu gwiazda:\n";
            std::cout << "Centralny wierzchołek: " << wierzcholek.podciag << " (ID: " << wierzcholek.idSekwencji << ", Pozycja: " << wierzcholek.pozycja << ")\n";
            std::cout << "Sąsiedzi:\n";
            for (const auto& sasiad : sasiedzi) {
                std::cout << "  Podciąg: " << sasiad.podciag << " (ID: " << sasiad.idSekwencji << ", Pozycja: " << sasiad.pozycja << ")\n";
            }
            return;
        }
    }

    std::cout << "Nie znaleziono struktury typu gwiazda spełniającej warunki.\n";
}


int main() {

    auto sekwencje = wczytajFasta();
    auto jakosci = wczytajJakosc();

    int progJakosci;
    std::cout << "Podaj próg wiarygodności: ";
    while (!(std::cin >> progJakosci)) {
        std::cout << "Błędna wartość, spróbuj ponownie: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    int dlugoscPodciagu;
    std::cout << "Podaj długość podciągu (od 4 do 9): ";
    while (!(std::cin >> dlugoscPodciagu) || dlugoscPodciagu < 4 || dlugoscPodciagu > 9) {
        std::cout << "Błędna wartość, spróbuj ponownie: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    auto przetworzoneSekwencje = usunNiskiejJakosciNukleotydy(sekwencje, jakosci, progJakosci);
    auto graf = utworzGraf(przetworzoneSekwencje, dlugoscPodciagu);
    auto krawedzie = polaczWierzcholki(graf, dlugoscPodciagu);

    for (const auto& para : przetworzoneSekwencje) {
        const auto& id = para.first;
        const auto& przetworzona = para.second;

        std::cout << "ID: " << id << "\n";
        std::cout << "Sekwencja po przetworzeniu: " << przetworzona.sekwencja << "\n";
        std::cout << "Pozycje zachowanych nukleotydów: ";
        for (const auto& pozycja : przetworzona.pozycje) {
            std::cout << pozycja << " ";
        }
        std::cout << "\nPozycje i oceny jakości usuniętych nukleotydów: ";
        for (const auto& usuniety : przetworzona.usunieteNukleotydy) {
            std::cout << usuniety.first << "(" << usuniety.second << ") ";
        }
        std::cout << "\n\n";
    }
    wyswietlWierzcholkiGrafu(graf);

    for (const auto& krawedz : krawedzie) {
        std::cout << "Krawedz miedzy podciagami: "
                  << krawedz.wierzcholek1.podciag << " (ID: " << krawedz.wierzcholek1.idSekwencji << ", Pozycja: " << krawedz.wierzcholek1.pozycja << ")"
                  << " a "
                  << krawedz.wierzcholek2.podciag << " (ID: " << krawedz.wierzcholek2.idSekwencji << ", Pozycja: " << krawedz.wierzcholek2.pozycja << ")\n";
    }


    znajdzStruktureGwiazdy(graf, krawedzie);

    return 0;
}

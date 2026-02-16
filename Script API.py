from Bio import Entrez
import requests
import time
from datetime import datetime, timedelta
import os

# Configuration
Entrez.email = "robin.prigent@etu.univ-amu.fr"  # Obligatoire pour PubMed
MISTRAL_API_KEY = "ydv0RIWmA8mgHsyMveM86XLYGr7wGvb1"  # Remplace par ta clé Mistral
WEBHOOK_URL = "https://hook.eu1.make.com/xjs6uwedm6w3ufyu2j6becad80tfqcp0"  # Remplace par ton URL de webhook
LAST_RUN_FILE = "last_run.txt"  # Fichier pour sauvegarder la dernière exécution

def load_last_run_date():
    """Charge la date de la dernière exécution."""
    try:
        if os.path.exists(LAST_RUN_FILE):
            with open(LAST_RUN_FILE, "r") as f:
                last_run_date = f.read().strip()
                print(f"Date de la dernière exécution chargée depuis {LAST_RUN_FILE}: {last_run_date}")
                return last_run_date
        else:
            print(f"Fichier {LAST_RUN_FILE} non trouvé. Utilisation de la date par défaut (7 jours en arrière).")
    except Exception as e:
        print(f"Erreur lors de la lecture de {LAST_RUN_FILE}: {e}. Utilisation de la date par défaut.")
    return (datetime.now() - timedelta(days=7)).strftime("%Y/%m/%d")

def save_last_run_date():
    """Sauvegarde la date actuelle comme dernière exécution."""
    try:
        with open(LAST_RUN_FILE, "w") as f:
            current_date = datetime.now().strftime("%Y/%m/%d")
            f.write(current_date)
        print(f"Date de la dernière exécution sauvegardée dans {LAST_RUN_FILE}: {current_date}")
    except Exception as e:
        print(f"Erreur lors de l'écriture dans {LAST_RUN_FILE}: {e}")

def search_pubmed(query, retmax=20):
    """Recherche les nouveaux articles depuis la dernière exécution."""
    last_run = load_last_run_date()
    date_today = datetime.now().strftime("%Y/%m/%d")
    print(f"Recherche des articles publiés entre {last_run} et {date_today}")
    search_query = f'{query} AND ("{last_run}"[Date - Publication] : "{date_today}"[Date - Publication]) NOT "review"[Publication Type]'
    handle = Entrez.esearch(db="pubmed", term=search_query, retmax=retmax, sort="pub_date")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_abstract(article_id):
    """Récupère l'abstract d'un article donné."""
    handle = Entrez.efetch(db="pubmed", id=article_id, rettype="abstract", retmode="text")
    abstract = handle.read()
    handle.close()
    return abstract

def analyze_with_mistral(text):
    headers = {
        "Authorization": f"Bearer {MISTRAL_API_KEY}",
        "Content-Type": "application/json"
    }
    payload = {
        "model": "mistral-tiny",
        "messages": [{
            "role": "user",
            "content": f"""
            Analyse cet abstract en français et extrais :
            1) Objectif de l'étude.
            2) Méthodes (type de LNPs, modèle animal, etc.).
            3) Résultats principaux.
            4) Innovations ou limites.
            Abstract: {text}
            """
        }]
    }
    response = requests.post("https://api.mistral.ai/v1/chat/completions", json=payload, headers=headers)
    response_data = response.json()

    if response.status_code != 200:
        print(f"Erreur {response.status_code}: {response.text}")
        return f"Erreur lors de l'appel à l'API: {response.text}"

    if "choices" not in response_data:
        print("Réponse inattendue de l'API :", response_data)
        return "Impossible d'extraire la réponse."

    return response_data["choices"][0]["message"]["content"]

def send_to_webhook(article_id, abstract, analysis):
    """Envoie les résultats à un webhook (ex: Make.com ou Zapier)."""
    data = {
        "article_id": article_id,
        "abstract": abstract,
        "analysis": analysis,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    response = requests.post(WEBHOOK_URL, json=data)
    print(f"Webhook response for article {article_id}: {response.status_code}")

def main():
    last_run_date = load_last_run_date()
    print(f"Date de la dernière exécution: {last_run_date}")

    query = '("chimeric antigen receptor T cells" OR "CAR-T" OR "CAR T") AND ("lipid nanoparticle" OR "LNP" OR "LPNs" OR "mRNA-LNPs")'
    article_ids = search_pubmed(query, retmax=20)

    if not article_ids:
        print("Aucun nouvel article trouvé.")
    else:
        print(f"Nouveaux articles trouvés : {len(article_ids)}")

        for article_id in article_ids:
            print(f"Traitement de l'article {article_id}...")
            abstract = fetch_abstract(article_id)
            analysis = analyze_with_mistral(abstract)
            send_to_webhook(article_id, abstract, analysis)
            time.sleep(10)  # Délai pour respecter les limites de l'API Mistral

    save_last_run_date()
    print("Analyse terminée. Résultats envoyés au webhook.")

if __name__ == "__main__":
    main()

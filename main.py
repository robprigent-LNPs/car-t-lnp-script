from Bio import Entrez
import requests
import time
from datetime import datetime, timedelta
import os

# Configuration
Entrez.email = "robin.prigent@etu.univ-amu.fr" 
MISTRAL_API_KEY = os.getenv("MISTRAL_API_KEY") 
WEBHOOK_URL = os.getenv("WEBHOOK_URL") 
LAST_RUN_FILE = "last_run.txt" 
PROCESSED_IDS_FILE = "processed_ids.txt" # NOUVEAU : Fichier pour stocker les IDs

def load_last_run_date():
    """Charge la date de la dernière exécution."""
    try:
        if os.path.exists(LAST_RUN_FILE):
            with open(LAST_RUN_FILE, "r") as f:
                return f.read().strip()
    except Exception:
        pass
    return (datetime.now() - timedelta(days=7)).strftime("%Y/%m/%d")

def save_last_run_date():
    """Sauvegarde la date actuelle."""
    try:
        with open(LAST_RUN_FILE, "w") as f:
            f.write(datetime.now().strftime("%Y/%m/%d"))
    except Exception as e:
        print(f"Erreur sauvegarde date: {e}")

# --- NOUVELLES FONCTIONS POUR LES IDs ---
def load_processed_ids():
    """Charge les IDs déjà traités pour éviter les doublons."""
    if os.path.exists(PROCESSED_IDS_FILE):
        with open(PROCESSED_IDS_FILE, "r") as f:
            return set(line.strip() for line in f)
    return set()

def save_processed_id(article_id):
    """Ajoute un ID traité au fichier."""
    try:
        with open(PROCESSED_IDS_FILE, "a") as f:
            f.write(f"{article_id}\n")
    except Exception as e:
        print(f"Erreur sauvegarde ID: {e}")
# ----------------------------------------

def search_pubmed(query, retmax=20):
    last_run = load_last_run_date()
    date_today = datetime.now().strftime("%Y/%m/%d")
    
    # On cherche large : depuis la dernière exécution jusqu'à aujourd'hui
    print(f"Recherche entre {last_run} et {date_today}")
    
    # Utilisation de [Date - Entrez] souvent plus précis pour les nouveautés que [Date - Publication]
    search_query = f'{query} AND ("{last_run}"[Date - Entrez] : "{date_today}"[Date - Entrez]) NOT "review"[Publication Type]'
    
    handle = Entrez.esearch(db="pubmed", term=search_query, retmax=retmax, sort="pub_date")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_abstract(article_id):
    handle = Entrez.efetch(db="pubmed", id=article_id, rettype="abstract", retmode="text")
    abstract = handle.read()
    handle.close()
    return abstract

def analyze_with_mistral(text):
    # (Ton code inchangé ici...)
    headers = {
        "Authorization": f"Bearer {MISTRAL_API_KEY}",
        "Content-Type": "application/json"
    }
    payload = {
        "model": "mistral-tiny",
        "messages": [{
            "role": "user",
            "content": f"Analyse cet abstract en français (bref): 1) Objectif 2) Méthodes 3) Résultats. Abstract: {text}"
        }]
    }
    try:
        response = requests.post("https://api.mistral.ai/v1/chat/completions", json=payload, headers=headers)
        if response.status_code == 200:
            return response.json()["choices"][0]["message"]["content"]
    except Exception as e:
        return f"Erreur API: {e}"
    return "Erreur analyse"

def send_to_webhook(article_id, abstract, analysis):
    # (Ton code inchangé ici...)
    data = {
        "article_id": article_id,
        "analysis": analysis,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    try:
        requests.post(WEBHOOK_URL, json=data)
    except Exception as e:
        print(f"Erreur webhook: {e}")

def main():
    last_run_date = load_last_run_date()
    processed_ids = load_processed_ids() # Chargement des anciens IDs
    
    print(f"Dernière exécution : {last_run_date}")
    
    query = '("chimeric antigen receptor T cells" OR "CAR-T" OR "CAR T") AND ("lipid nanoparticle" OR "LNP" OR "mRNA-LNPs")'
    article_ids = search_pubmed(query, retmax=20)

    if not article_ids:
        print("Aucun article trouvé dans la plage de dates.")
    else:
        new_articles_count = 0
        for article_id in article_ids:
            # VÉRIFICATION DOUBLON
            if article_id in processed_ids:
                print(f"Article {article_id} déjà traité. Ignoré.")
                continue
            
            print(f"Traitement du nouvel article {article_id}...")
            abstract = fetch_abstract(article_id)
            
            # Petite sécurité si pas d'abstract
            if not abstract:
                print("Pas d'abstract disponible.")
                continue

            analysis = analyze_with_mistral(abstract)
            send_to_webhook(article_id, abstract, analysis)
            
            # Sauvegarde immédiate de l'ID pour ne pas le refaire si le script plante après
            save_processed_id(article_id)
            new_articles_count += 1
            
            time.sleep(5) 

        print(f"Analyse terminée. {new_articles_count} nouveaux articles envoyés.")

    save_last_run_date()

if __name__ == "__main__":
    main()

#!/usr/bin/env python3

"""The equivalent of 'gcloud auth application-default print-access-token' for a service
account with JSON credentials, without requiring the Google Cloud SDK."""

import os
import google.auth
import google.auth.transport.requests
from google.oauth2 import service_account

credentials_json = os.getenv('GOOGLE_APPLICATION_CREDENTIALS', '/gsa-key/key.json')
credentials = service_account.Credentials.from_service_account_file(
    credentials_json, scopes=['https://www.googleapis.com/auth/cloud-platform']
)
credentials.refresh(google.auth.transport.requests.Request())
print(credentials.token)

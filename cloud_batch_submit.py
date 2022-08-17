#!/usr/bin/env python3

import argparse
import subprocess
import textwrap
import uuid

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # See Cloud Build variable substitutions.
    parser.add_argument("--location", required=True)
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--tag-name", required=True)
    # See cuking.cu for argument help.
    parser.add_argument("--input-uri", required=True)
    parser.add_argument("--output-uri", required=True)
    parser.add_argument("--requester-pays-project", required=True)
    parser.add_argument("--kin-threshold", type=float, required=True)
    parser.add_argument("--split-factor", type=int, required=True)

    args = parser.parse_args()

    batch_job_json = textwrap.dedent(
        """\
        {{
          "taskGroups": [
            {{
              "taskSpec": {{
                "runnables": [
                  {{
                    "script": {{
                      "text": "sudo docker run --name cuking --gpus all {location}-docker.pkg.dev/{project_id}/images/cuking:{tag_name} cuking --input_uri={input_uri} --output_uri={output_uri} --requester_pays_project={requester_pays_project} --kin_threshold={kin_threshold} --split_factor={split_factor} --shard_index=${{BATCH_TASK_INDEX}}"
                    }}
                  }}
                ],
                "computeResource": {{
                  "cpuMilli": 12000,
                  "memoryMib": 87040
                }},
                "maxRunDuration": "36000s"
              }},
              "taskCount": {task_count}
            }}
          ],
          "allocationPolicy": {{
            "instances": [
              {{
                "instanceTemplate": "cuking-instance-template"
              }}
            ]
          }},
          "logsPolicy": {{
            "destination": "CLOUD_LOGGING"
          }}
        }}
        """
    ).format(**vars(args), task_count=args.split_factor * (args.split_factor + 1) // 2)

    JSON_FILENAME = 'batch_job.json'
    with open(JSON_FILENAME, 'w') as f:
        print(batch_job_json, file=f)

    # Use a UUID in the job name to avoid collisions.
    job_name = f'cuking-{uuid.uuid4()}'
    cmd = [
        'gcloud',
        'beta',
        'batch',
        'jobs',
        'submit',
        job_name,
        f'--location={args.location}',
        f'--config={JSON_FILENAME}',
    ]
    print(f'Submitting job:\n    {" ".join(cmd)}')
    subprocess.run(cmd, check=True)

    print(
        f'\nTo check the status of the job, run:'
        f'\n    gcloud beta batch jobs describe {job_name} --location={args.location}'
    )
